#= Implements optical scattering under the reduced Rayleigh equations
# Effectively solves Maxwell's equations for a singularly polarized wave
# on a partially symmetric rough surface. Satisfying the boundary conditions

# The reduced Rayleigh equations are a set of coupled integral equations
# which assume that far field scattering conditions (singular direction, up/down in 1D)
# can be used all the way down to the rough surface boundary even though for strongly
# rough surfaces one can get multiple scattering events.
=#


using FFTW
using Plots
using BenchmarkTools
using LinearAlgebra




@enum Polarization p s
struct RayleighParams

    FT::FFTW.cFFTWPlan{ComplexF64,-1,true,1,Tuple{Int64}} # Planned Fourier transform of surface points
    M_pre::Array{ComplexF64,3} # Precomputed Mpq coefficients
    N_pre::Matrix{ComplexF64} # Precomputed Npk coefficients

    xs::Vector{Float64}  # Surface points
    ps::Vector{Float64}  # Scattered wave numbers
    qs::Vector{Float64}  # Scattered wave numbers
    wq::Vector{Float64}  # Weights of the quadrature in q_n

    ν::Polarization # Polarization [p, s]
    ε::ComplexF64 # Permittivity of the scattering medium
    μ::ComplexF64 # Permeability of the scattering medium
    λ::Float64   # wavelength in [μm]
    ω::Float64   # Angular frequency
    k::Float64   # Parallel component wave number
    ki::Int      # Looked up index of the incoming wave number

    # Sizings
    Q::Float64   # Truncated wave number, q = (-∞,∞) -> q_n = (-Q_mult / 2, Q_mult / 2)
    Δq::Float64  # Wave number spacing [1/m]
    Lx::Float64  # Surface length [m]
    Δx::Float64  # Surface point spacing [m]

    Ni::Int      # Order of surface power expansion

    # Constructor
    RayleighParams(; ν::Polarization=p, ε=2.25, μ=1.0, λ=600e-9, Q_mult=4, Nq=127, L=10.0e-6, Ni=10, θ0=45) = begin
        #=
        All lengths are in units of μm
        =#
        display("Initialising variables...")

        c0 = 299_792_458.0 # [m/s], Speed of light in vacuum
        # c0 = 1.0 # [m/s], Speed of light in vacuum, natural units

        K = 2π / λ
        ω = c0 * K
        k = K * sind(θ0) # Parallel component wave number

        # All variables scaled such that ω/c0 = 1

        # Assertions and warnings
        @assert L / λ > 10.0 "Surface length must be much larger than the wavelength, L ≫ λ, but is L:$L and λ:$λ."
        @assert Q_mult > 2 "Q_mult must be greater than 2, but is $Q_mult. ¤ is recommended."
        @assert Nq > 2 "Nq must be greater than 2, but is $Nq."

        Nx = 2 * Nq

        # Q = Q_mult * ω / c0 # Truncated wave number
        Q = Q_mult # Truncated wave number, divided by ω / c0
        Δq = Q / Nq

        ps = -Q/2:Δq:Q/2
        qs = Q/2:-Δq:-Q/2

        # Find the nearest index of the incoming wave number
        ki = searchsortedfirst(qs, k, rev=true)
        k_new = qs[ki]

        wq = ones(size(qs))
        wq[1] = wq[end] = 3.0 / 8.0
        wq[2] = wq[end-1] = 7.0 / 6.0
        wq[3] = wq[end-2] = 23.0 / 24.0

        Lx = L * ω / c0 # Surface length, scaled up by ω / c0, since reciprocal space is scaled down by ω / c0
        Δx = Lx / Nx
        xs = -Lx/2:Δx:Lx/2

        ε_new = ε + 1e-6im # Add a small imaginary part to avoid singularities

        ## Precalculate invariant coefficients:
        display("Preconditioning solver")
        N_pre = Matrix{ComplexF64}(undef, length(ps), Ni + 1)
        M_pre = Array{ComplexF64,3}(undef, length(ps), length(qs), Ni + 1)

        α(q::Float64)::ComplexF64 = √complex(ε * μ - q^2)

        # Assumed μ0 = ε0 = 1
        α0(q::Float64)::ComplexF64 = √complex(1.0 - q^2)

        κ = ν == p ? ε : μ

        display("Precalculating invariant coefficients in M.")
        Mpq_invariant(p::Float64, q::Float64, n::Int)::ComplexF64 = (
            (-1.0im)^n / factorial(n) * (
                (p + κ * q) * (p - q) * (α(p) - α0(q))^(n - 1) +
                (α(p) + κ * α0(q)) * (α(p) - α0(q))^n)
        )
        for I in CartesianIndices(M_pre)
            i, j, n = Tuple(I)
            M_pre[i, j, n] = Mpq_invariant(ps[i], qs[j], n - 1)
        end

        display("Precalculating invariant coefficients in N.")
        Npk_invariant(p::Float64, n::Int)::ComplexF64 = (
            (-1.0im)^n / factorial(n) * (
                (p + κ * k) * (p - k) * (α(p) + α0(k))^(n - 1) +
                (α(p) - κ * α0(k)) * (α(p) + α0(k))^n)
        )

        for I in CartesianIndices(N_pre)
            i, n = Tuple(I)
            N_pre[i, n] = Npk_invariant(ps[i], n - 1)
        end

        new(plan_fft!(xs), M_pre, N_pre,
            xs, ps, qs, wq,
            ν, ε_new, μ, λ, ω, k_new, ki,
            Q, Δq, Lx, Δx, Ni)
    end
end

function show_params(rp::RayleighParams)
    names = fieldnames(typeof(rp))
    for name in names
        field = getfield(rp, name)
        Base.print("$name:\t")
        display(field)
    end
end

get_θ(rp::RayleighParams) = asind(rp.k)

function gaussian_surface_gen(δ::Float64, a::Float64, rp::RayleighParams)
    W(x) = exp(-x^2 / a^2)
    g(k) = √π * a * exp(-(a * k / 2.0)^2)

    xs = rp.xs

    W_p = W.(xs) |> ComplexF64

    return rp.FT \ (
        (rp.FT * (randn(Float64, N) |> ComplexF64)) .*
        (P * W_p)) .* δ * 2π / sqrt(N)

end


@enum SurfType flat gauss singlebump
struct SurfPreAlloc
    # Preallocated steps in the calculations for a given surface
    Mpq::Matrix{ComplexF64}  # Matrix of the Mpq coefficients (A)
    Npk::Vector{ComplexF64}  # Vector of the Npk coefficients (b)
    Rqk::Vector{ComplexF64}  # Vector of the reflection solution (x = A \ b)

    Fys::Vector{ComplexF64}  # Fourier transform of surface heights, prealloc
    sFys::Vector{ComplexF64} # Shifted Fourier transform of surface heights, prealloc
    ys::Vector{Float64}  # Surface heights

    function SurfPreAlloc(rp::RayleighParams, surf_t::SurfType, δ::Float64, a::Float64)
        Lx = rp.Lx
        xs = rp.xs

        if surf_t == singlebump
            ys = δ * exp.(-xs .^ 2 / a^2)
        elseif surf_t == gauss

            # δ = RMS height of surface profile function
            # a = Autocorrelation length
            # TODO: Move code from SurfaceGen into here
        elseif surf_t == flat
            ys = zeros(rp.Nx + 1) # Flat surface
        end

        ys *= (ω / c0) # Scale surface heights by ω / c0

        ys = similar(xs)
        Fys = Vector{ComplexF64}(undef, rp.Nq + 1)
        sFys = similar(Fys)
        Mpq = zeros(ComplexF64, rp.Nq + 1, rp.Nq + 1)
        Npk = zeros(ComplexF64, rp.Nq + 1)
        Rqk = similar(Npk)

        new(Mpq, Npk, Rqk, Fys, sFys, ys)
    end
end

rp_test = RayleighParams(;
    ν=p,
    Nq=2^11,
    ε=2.25,
    L=10.0e-6,
    Q_mult=4,
    Ni=10
);



function solve(rp::RayleighParams, M_pre::Array{ComplexF64,3}, N_pre::Matrix{ComplexF64})

    for n in axes(M_pre, 3)
        display("n = $n")

        rp.Fys .= ys .^ n
        rp.FT * rp.Fys # In place FFT
        fftshift!(rp.sFys, rp.Fys) # No way to do shift in place apparently


        for I in CartesianIndices(rp.Mpq)
            i, j = Tuple(I)
            rp.Mpq[I] += M_pre[i, j, n] * rp.sFys[i+j-1]
        end

        for i in eachindex(rp.Npk)
            rp.Npk[i] -= N_pre[i, n] * rp.sFys[i+rp.ki-1]
        end
    end

    return rp.Mpq \ rp.Npk
end



function solve(rp::RayleighParams, reduction::Function=identity)
    # Solve the reduced Rayleigh equations for p-polarized light

    # Shorthand parameters
    κ = rp.ν == p ? rp.ε : rp.μ
    k = rp.k
    ki = rp.ki
    Ni = rp.Ni

    α(q::Float64)::ComplexF64 = √complex(rp.ε * rp.μ - q^2)

    # Assumed μ0 = ε0 = 1
    α0(q::Float64)::ComplexF64 = √complex(1.0 - q^2)

    Mpq_eval(p::Float64, q::Float64, n::Int)::ComplexF64 = (
        (-1.0im)^n / factorial(n) * (
            (p + κ * q) * (p - q) * (α(p) - α0(q))^(n - 1) +
            (α(p) + κ * α0(q)) * (α(p) - α0(q))^n)
    )
    Npk_eval(p::Float64, n::Int)::ComplexF64 = (
        (-1.0im)^n / factorial(n) * (
            (p + κ * k) * (p - k) * (α(p) + α0(k))^(n - 1) +
            (α(p) - κ * α0(k)) * (α(p) + α0(k))^n)
    )

    for n in 0:Ni
        display("n = $n")

        rp.Fys .= rp.ys .^ n
        rp.FT * rp.Fys # In place FFT
        fftshift!(rp.sFys, rp.Fys) # No way to do shift in place apparently

        for I in CartesianIndices(rp.Mpq)
            i, j = Tuple(I)
            rp.Mpq[I] += Mpq_eval(rp.ps[i], rp.qs[j], n) * rp.sFys[i+j-1]
        end

        for i in eachindex(rp.Npk)
            rp.Npk[i] -= Npk_eval(rp.ps[i], n) * rp.sFys[i+ki-1]
        end
    end

    nans1 = findall(!isfinite, rp.Mpq)
    nans2 = findall(!isfinite, rp.Npk)

    if length(nans1) > 0 || length(nans2) > 0
        println("Mpq:")
        display(rp.Mpq[nans1])
        println("Mpq idxs:")
        display(nans1)

        println("Npk:")
        display(rp.Npk[nans2])
        println("Npk idxs:")
        display(nans2)

        display("Number of NaNs in Mpq: $(length(nans1))")
        display("Number of NaNs in Npk: $(length(nans2))")

        throw("Non-finite values in Mpq or Npk")
    end

    # Solve the system Δq / 2π ∑_q N⁺ₚ(p|q) * Rₚ(q|k) = N⁻ₚ(p|k) for R
    # rp.Rqk .= rp.Mpq \ rp.Npk

    return reduction(rp.Mpq \ rp.Npk)
end

mean_DRC(rp::RayleighParams, sol::Matrix{ComplexF64}) = begin
    # Calculate the mean differential reflection coefficient
    # Factor of ω / c0 removed due to scaling

    ω = rp.ω
    c0 = rp.c0
    qs = rp.qs[rp.qs.<ω/c0.&&rp.qs.>-ω/c0]
    ks = rp.ks

    retval = Matrix{Float64}(undef, length(qs), length(ks))

    for (i, q) in enumerate(qs)
        for (j, k) in enumerate(ks)
            retval[i, j] = 1.0 / rp.Lx / 2π * cos(θ(q, ω, c0))^2 / cos(θ(k, ω, c0)) * abs2(sol[i, j])
        end
    end
    return retval
end

Nq = 2^11

silver = -17.5 + 0.48im
glass = 2.25

rp_p = RayleighParams(;
    ν=p,
    Nq=Nq,
    ε=glass,
    L=10.0e-6,
    Q_mult=4,
    Ni=10
);


reduction = x::Vector{ComplexF64} -> maximum(abs.(x))^2

# show_params(rp_p)
@time sol_p = solve(rp_p, reduction);

sol_p


# sum(abs2, sol_p)
rp_s = RayleighParams(;
    ν=s,
    Nq=Nq,
    ε=glass,
    L=10.0e-6,
    Q_mult=4,
    Ni=10
);
# show_params(rp_s)
@time sol_s = solve(rp_s);
θs_new_p = asind.(rp_p.ks)
θs_new_s = asind.(rp_s.ks)

refl_p = (maximum(abs.(sol_p), dims=1) |> vec) .^ 2
refl_s = (maximum(abs.(sol_s), dims=1) |> vec) .^ 2

display(refl_s[1])



ε = glass

r_analytical = abs2.((1.0 .- .√(ε .- rp_p.ks .^ 2)) ./ (1 .+ .√(ε .- rp_p.ks .^ 2)))

brewster = 56

bi = findfirst(refl_p .== minimum(refl_p))
display("Brewster angle $brewster vs. $(θs_new_p[bi])")

plot(θs_new_p, refl_p, label="ν = p", marker=(:circle, 2));
plot!(θs_new_s, refl_s, label="ν = s", marker=(:circle, 2));
plot!(θs_new_p, r_analytical, label="Analytical Fresnel")
plot!([0.0, 90.0], [1.0, 1.0], linestyle=:dash, linecolor=:black, linewidth=1, label=nothing);
plot!([90.0, 90.0], [0.0, 1.0], linestyle=:dash, linecolor=:black, linewidth=1, label=nothing);


xlims!(0, 90);
ylims!(0, 1.00);
xlabel!("Θ [°]");
ylabel!("Fresnel reflection coefficient")