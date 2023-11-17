module RayleighSolver

using RayleighSetup
using FFTW

export α, α0, solve, pre_M_invariant!, pre_N_invariant!, solve_pre!

function α(q::Float64, εμ::ComplexF64)::ComplexF64
    # Calculate the α parameter for a given p
    return sqrt(εμ - q^2)
end
function α0(q::Float64)::ComplexF64
    # Calculate the α0 parameter for a given q
    return sqrt(complex(1.0 - q^2))
end
function M_ker(p::Float64, q::Float64, κ::ComplexF64, α::ComplexF64, α0::ComplexF64, n::Int)::ComplexF64
    # Calculate the kernel of the Mpq matrix
    Δα = α - α0
    return (-1.0im)^n / factorial(n) * (
        (p + κ * q) * (p - q) * Δα^(n - 1) +
        (α + κ * α0) * Δα^n
    )
end
function N_ker(p::Float64, k::Float64, κ::ComplexF64, α::ComplexF64, α0::ComplexF64, n::Int)::ComplexF64
    Δα = α + α0
    return (-1.0im)^n / factorial(n) * (
        (p + κ * k) * (p - k) * Δα^(n - 1) +
        (α - κ * α0) * Δα^n
    )
end
function M_invariant!(M::AbstractMatrix{ComplexF64}, ps::Vector{Float64}, qs::Vector{Float64}, κ::ComplexF64, εμ::ComplexF64, n::Int)::Nothing

    @inbounds for (i, p) in enumerate(ps)
        for (j, q) in enumerate(qs)
            M[i, j] = M_ker(p, q, κ, α(p, εμ), α0(q), n)
        end
    end
    return nothing
end
function N_invariant!(N::AbstractVector{ComplexF64}, ps::Vector{Float64}, k::Float64, κ::ComplexF64, εμ::ComplexF64, n::Int)::Nothing

    @inbounds for (i, p) in enumerate(ps)
        N[i] = N_ker(p, k, κ, α(p, εμ), α0(k), n)
    end
    return nothing
end
function M_invariant(ps::Vector{Float64}, qs::Vector{Float64}, κ::ComplexF64, εμ::ComplexF64, n::Int)::Matrix{ComplexF64}
    M = Matrix{ComplexF64}(undef, length(ps), length(qs))
    M_invariant!(M, ps, qs, κ, εμ, n)
    return M
end


function N_invariant(ps::Vector{Float64}, k::Float64, εμ::ComplexF64, κ::ComplexF64, n::Int)::Vector{ComplexF64}
    N = similar(ps, ComplexF64)
    N_invariant!(N, ps, k, εμ, κ, n)
    return N
end

function pre_M_invariant!(M::Array{ComplexF64,3}, rp::RayleighParams)::Nothing
    # Calculate the surface invariant part of the Mpq matrix
    # Invariant under surface change), but NOT incident angle θ0

    ps = rp.ps
    qs = rp.qs
    κ = rp.ν == p ? rp.ε : rp.μ
    εμ = rp.ε * rp.μ

    @inbounds for n in axes(M, 3)
        Mn = view(M, :, :, n)
        M_invariant!(Mn, ps, qs, κ, εμ, n - 1)
    end
    return nothing
end

function pre_N_invariant!(N::Matrix{ComplexF64}, rp::RayleighParams, k::Float64)::Nothing
    # Calculate the surface invariant part of the Npk vector
    # Invariant under surface change, but NOT incident angle θ0
    # Input k must be the nearest value in qs for correct results
    ps = rp.ps
    κ = rp.ν == p ? rp.ε : rp.μ
    εμ = rp.ε * rp.μ

    @inbounds for n in axes(N, 2)
        Nn = view(N, :, n)
        N_invariant!(Nn, ps, k, εμ, κ, n - 1)
    end
    return nothing
end

function solve_pre!(sp::SurfPreAlloc, rp::RayleighParams,
    M_pre::Array{ComplexF64,3}, N_pre::Matrix{ComplexF64},
    ki::Int)::Nothing
    sp.Mpq .= 0.0
    sp.Npk .= 0.0

    for n in axes(M_pre, 3)
        sp.Fys .= sp.ys .^ (n - 1)
        rp.FT * sp.Fys # In place FFT
        fftshift!(sp.sFys, sp.Fys)

        @inbounds for I in eachindex(IndexCartesian(), sp.Mpq)
            i, j = Tuple(I)
            sp.Mpq[i, j] += M_pre[i, j, n] * sp.sFys[i+j-1]
        end

        @inbounds for i in eachindex(sp.Npk)
            sp.Npk[i] -= N_pre[i, n] * sp.sFys[i+ki-1]
        end
    end

    sp.R .= sp.Mpq \ sp.Npk
    return nothing
end

function solve(rp::RayleighParams, ys::Vector{Float64}, ki::Int)::Vector{ComplexF64}
    # Allocating version of the solver
    # Solve for the reflection coefficient R
    #   version of solver without pre-calculated M and N
    ps = rp.ps
    qs = rp.qs
    κ = rp.ν == p ? rp.ε : rp.μ
    εμ = rp.ε * rp.μ
    Ni = rp.Ni
    k = qs[ki]

    Mpq = zeros(ComplexF64, length(ps), length(qs))  # Matrix of the Mpq coefficients (A)
    Npk = zeros(ComplexF64, length(ps))  # Vector of the Npk coefficients (b)
    Fys = similar(ys, ComplexF64)  # Fourier transform of surface heights, prealloc

    for n in 0:Ni
        Fys .= ys .^ n
        rp.FT * Fys # In place FFT
        Fys .= fftshift(Fys)

        for (i, p) in enumerate(ps)
            for (j, q) in enumerate(qs)
                Mpq[i, j] += M_ker(p, q, κ, α(p, εμ), α0(q), n) * Fys[i+j-1]
            end
            Npk[i] -= N_ker(p, k, κ, α(p, εμ), α0(k), n) * Fys[i+ki-1]
        end
    end

    return Mpq \ Npk
end
end # module RayleighSolver
