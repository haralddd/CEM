
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

function M_invariant!(M::Array{ComplexF64,3}, rp::RayleighParams)::Nothing
    # Calculate the surface invariant part of the Mpq matrix
    # Invariant under surface change), but NOT incident angle θ0

    ps = rp.ps
    qs = rp.qs
    κ = rp.ν == p ? rp.ε : rp.μ
    εμ = rp.ε * rp.μ

    for n in axes(M, 3)
        for i in eachindex(ps)
            p = ps[i]
            for j in eachindex(qs)
                q = qs[j]
                M[i, j, n] = M_ker(p, q, κ, α(p, εμ), α0(q), n - 1)
            end
        end
    end
    return nothing
end

function N_invariant!(N::Matrix{ComplexF64}, rp::RayleighParams, k::Float64)::Nothing
    # Calculate the surface invariant part of the Npk vector
    # Invariant under surface change, but NOT incident angle θ0
    # Input k must be the nearest value in qs for correct results
    ps = rp.ps
    κ = rp.ν == p ? rp.ε : rp.μ
    εμ = rp.ε * rp.μ

    for n in axes(N, 2)
        for i in eachindex(ps)
            p = ps[i]
            N[i, n] = N_ker(p, k, κ, α(p, εμ), α0(k), n - 1)
        end
    end
    return nothing
end

function solve!(sp::SurfPreAlloc, rp::RayleighParams,
    M_pre::Array{ComplexF64,3}, N_pre::Matrix{ComplexF64},
    ki::Int)::Nothing

    # Calculates the preallocated surface integral
    # sp is the preallocated surface struct, containing the preallocated output
    # rp is the RayleighParams struct, containing the constant parameters for the calculation
    # Matrix M_pre is preallocated and contains the invariant part of the Mpq matrix
    # Matrix N_pre is preallocated and contains the invariant part of the Npk vector
    # ki is the index of the nearest value in qs to the incident angle θ0

    # Stores the result in sp.Npk

    sp.Mpq .= 0.0
    sp.Npk .= 0.0

    for n in axes(M_pre, 3)
        sp.Fys .= sp.ys .^ (n - 1)
        rp.FT * sp.Fys # In place FFT
        fftshift!(sp.sFys, sp.Fys)

        for i in axes(sp.Mpq, 1), j in axes(sp.Mpq, 2)
            sp.Mpq[i, j] += M_pre[i, j, n] * sp.sFys[i+j-1]
        end

        for i in eachindex(sp.Npk)
            sp.Npk[i] -= N_pre[i, n] * sp.sFys[i+ki-1]
        end
    end

    sp.Npk .= sp.Mpq \ sp.Npk
    return nothing
end