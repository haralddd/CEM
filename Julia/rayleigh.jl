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
using StaticArrays
using LinearAlgebra
# using Statistics


@enum Polarization p s
@enum SurfType flat gauss

v = rand(1024);
F = plan_rfft(v)
typeof(F)

vsz = SizedVector

struct RayleighParams{Nq,Ni}

    FT::FFTW.rFFTWPlan{Float64,-1,false,1,Tuple{Int64}} # Planned Fourier transform of surface points
    ζs::SizedVector{2 * Nq,Float64}  # Surface heights
    Fζs::SizedVector{Nq + 1,ComplexF64}  # Fourier transform of surface heights, prealloc
    ps::SizedVector{Nq,Float64}  # Scattered wave numbers
    qs::SizedVector{Nq,Float64}  # Scattered wave numbers


    ν::Polarization # Polarization [p, s]
    ε::Float64 # Permittivity of the scattering medium
    μ::Float64 # Permeability of the scaterring medium
    λ::Float64   # wavelength in [μm]
    ω::Float64   # Angular frequency
    k::Float64   # Parallel component wave number

    # Sizings
    Q::Float64   # Truncated wave number, q = (-∞,∞) -> q_n = (-Q_mult / 2, Q_mult / 2)
    Δq::Float64  # Wave number spacing [1/m]
    # wq::Float64  # Weights of the quadrature in q_n

    Lx::Float64  # Surface length [m]
    Δx::Float64  # Surface point spacing [m]

    RayleighParams(; ν::Polarization=p, surf_t::SurfType=flat, ε=2.25, μ=1.0, λ=600e-9, Q_mult=4, Nq=127, L=10.0e-6, Ni=10, θ0=45) = begin
        #=
        All lengths are in units of μm
        =#

        c0 = 299_792_458.0 # [m/s], Speed of light in vacuum
        # c0 = 1.0 # [m/s], Speed of light in vacuum, natural units

        K = 2π / λ
        ω = c0 * K

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

        # ks = 0.0:Δq:1.0
        # ks = sind.(0.0:0.5:90.0)

        kis = [searchsortedfirst(qs, ks[i], rev=true) for i in eachindex(ks)] |> unique
        display(kis)
        ks_new = qs[kis]
        # wq = ones(size(qs))
        # wq[1] = wq[end] = 3.0 / 8.0
        # wq[2] = wq[end-1] = 7.0 / 6.0
        # wq[3] = wq[end-2] = 23.0 / 24.0

        Lx = L * ω / c0 # Surface length, scaled up by ω / c0, since reciprocal space is scaled down by ω / c0
        Δx = Lx / Nx
        xs = -Lx/2:Δx:Lx/2
        ζs = SizedVector{Nx}(Vector{ComplexF64}(undef, Nx))

        if false #surf_t == gauss
            ζs = cos.(4 * 2π .* xs ./ Lx) * 1e-9
        elseif surf_t == flat
            ζs = zeros(ComplexF64, Nx) # Flat surface
        end

        ζs *= (ω / c0) # Scale surface heights by ω / c0

        new{Nq,Ni}(plan_fft(ζs), ζs, ps, qs,
            ν, ε, μ, λ, ω, K,
            Q, Δq, wq, Lx, Δx)
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

function solve(rp::RayleighParams)
    # Solve the reduced Rayleigh equations for p-polarized light

    ε = rp.ε + 1e-6im # Permittivity of the medium, added small absorption to avoid singularities
    μ = rp.μ # Permeability of the medium

    κ = rp.ν == p ? ε : μ # Polarization dependent factor

    ω = rp.ω # Angular frequency
    c0 = rp.c0 # Speed of light in vacuum

    ps = rp.ps # Scattered ∥ wave numbers
    qs = rp.qs # Scattered ∥ wave numbers

    k = rp.k # Incoming ∥ wave number
    ki = rp.ki
    ζs = rp.ζs # Surface heights

    FT = rp.FT

    Nq = length(qs)
    Np = length(ps)
    Nk = length(ks)
    Nx = length(ζs)
    Ni = rp.Ni # Order of surface power expansion

    Mpq = zeros(ComplexF64, Np, Nq) # M⁺(p|q) matrix
    Rqk = Matrix{ComplexF64}(undef, Nq, Nk) # Vector of solution R(q|k) vectors
    Mpk = zeros(ComplexF64, Np, Nk) # Vector of M⁻(p|k) RHS vectors


    α(q::Float64)::ComplexF64 = √complex(ε * μ - q^2)

    # Assumed μ0 = ε0 = 1
    α0(q::Float64)::ComplexF64 =
        abs(q) > 1.0 ?
        1.0im * √(q^2 - 1.0) :
        √complex(1.0 - q^2)

    # Symmetrize gs access -NQ,...,-1, 0, 1,...,NQ, by p-q indexing 
    function pq_idx(i::Int, j::Int)::Int
        idx = -g_size + i + j - 1
        return (idx < 0) ? -idx + 1 : idx + 1
    end

    display("Evaluate M⁺(p,q) matrix: ")
    for m in 0:Ni
        mul!(FTζ, rp.FT_plan, ζs)
        for (i, p) in enumerate(ps)
            for (j, q) in enumerate(qs)
                denom = α(p) - α0(q)
                Mpq[i, j] += (-1.0im)^n * gs[pq_idx(i, j), m] / factorial(n) * (
                    (p + κ * q) * (p - q) * denom^(n - 1) +
                    (α(p) + κ * α0(q)) * denom^n
                )
            end

            denom = α(p) + α0(k)

            Mpk[i] -= (-1.0im)^n * FTζ[pq_idx(i, ki)] / factorial(n) * (
                (p + κ * k) * (p - k) * denom^(n - 1) +
                (α(p) - κ * α0(k)) * denom^n
            )
        end
    end

    display("Evaluate M⁻ₚ(p,k) vector: ")

    @time for m in axes(gs, 2)
        n = m - 1

    end

    nans1 = findall(!isfinite, Mpq)
    nans2 = findall(!isfinite, Mpk)

    if length(nans1) > 0 || length(nans2) > 0
        println("Mpq:")
        display(Mpq[nans1])
        println("Mpq idxs:")
        display(nans1)

        println("Mpk:")
        display(Mpk[nans2])
        println("Mpk idxs:")
        display(nans2)

        display("Number of NaNs in Mpq: $(length(nans1))")
        display("Number of NaNs in Mpk: $(length(nans2))")

        throw("Non-finite values in Mpq or Mpk")
    end

    # Solve the system Δq / 2π ∑_q N⁺ₚ(p|q) * Rₚ(q|k) = N⁻ₚ(p|k) for R
    display("Solve Ax = b: ")
    @time for i in eachindex(ks)
        Rqk[:, i] = Mpq \ Mpk[:, i]
    end

    return Rqk
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

Nq = 2^10
display(Nq)

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



# show_params(rp_p)
@time sol_p = solve(rp_p);



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