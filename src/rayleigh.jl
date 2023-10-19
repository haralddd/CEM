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

@enum Polarization p s
@enum SurfType flat gauss

function is_power_two(n::Int)::Bool
    # Check if n is a power of 2
    (n & (n - 1)) == 0
end

signal = sin.(range(0, stop=10π, length=1000))
plot(signal)
fft(signal) |> fftshift .|> abs |> plot
rfft(signal) |> fftshift .|> abs |> plot

struct RayleighParams
    ν::Polarization # Polarization [p, s]
    c0 # [m/s], Speed of light in vacuum
    κ0  # ν: ϵ₀ / μ₀
    κ::ComplexF64   # Polarization dependent constant
    λ   # Angular frequency

    # Sizings
    Q   # Truncated wave number, q = (-∞,∞) -> q_n = (-Q_mult / 2, Q_mult / 2)
    Nq  # Number of wave numbers
    Δq  # Wave number spacing [1/m]
    wq  # Weights of the quadrature in q_n
    # Qis  # Integer array from p - q to Qi[i, j] in the Fourier transform array

    Nx  # Number of surface points
    Lx  # Surface length [m]
    Δx  # Surface point spacing [m]

    Ni  # Order of surface power expansion

    ps   # Scattered wave numbers
    qs   # Scattered wave numbers
    ks   # Incoming wave numbers
    ζs   # Surface heights

    FT_plan # Planned Fourier transform of surface points

    RayleighParams(; ν::Polarization=p, surf_t::SurfType=flat, κ0=1.0, κ=2.25, λ=400e-9, Q_mult=4, Nq=128, L=1.0e-4, Ni=10) = begin
        c0 = 299_792_458.0
        # λ = c0 / 2πω
        K = 2π / λ
        ω = c0 * K

        # Scale all variables such that ω/c0 = 1

        # Assertions and warnings
        @assert L / λ > 100.0 "Surface length must be much larger than 1 wavelength, Lx ≫ λ, but is Lx:$Lx and λ:$λ."
        @assert Q_mult > 2 "Q_mult must be greater than 2, but is $Q_mult."
        @assert Nq > 2 "Nq must be greater than 2, but is $Nq."



        Nx = 2 * Nq

        Q = Q_mult # Truncated wave number
        Δq = Q / (Nq - 1.0)

        # Non-inclusive, such that Nq splits the range evenly.
        ps = -Q/2:Δq:Q/2
        qs = Q/2:-Δq:-Q/2
        ks = -1.0:Δq:1.0

        wq = ones(Nq + 1)
        wq[1] = wq[end] = 3.0 / 8.0
        wq[2] = wq[end-1] = 7.0 / 6.0
        wq[3] = wq[end-2] = 23.0 / 24.0

        Lx = L * (ω / c0)
        Δx = Lx / (Nx - 2.0)
        xs = -Lx/2:Δx:Lx/2

        if false #surf_t == gauss
            ζs = cos.(4 * 2π .* xs ./ Lx) * 1e-9
        elseif surf_t == flat
            ζs = zeros(Nx - 2) # Flat surface
        end
        # Rescale ζ to be in units of λ
        ζs *= (ω / c0)

        # Qis = Matrix{Int}(undef, Nq + 1, Nq + 1)

        new(ν, c0, κ0, κ, λ, Q, Nq, Δq, wq,
            Nx, Lx, Δx, Ni, ps, qs, ks, ζs,
            plan_rfft(ζs))
    end
end

function print(rp::RayleighParams)
    names = fieldnames(typeof(rp))
    for name in names
        field = getfield(rp, name)
        Base.print("$name:\t")
        display(field)
    end
end

function solve_p(rp::RayleighParams)
    # Solve the reduced Rayleigh equations for p-polarized light

    ε0 = rp.κ0 # Permittivity of free space
    ε = rp.κ # Permittivity of the medium

    ps = rp.ps # Scattered wave numbers
    qs = rp.qs # Scattered wave numbers

    ks = rp.ks # Incoming wave numbers
    ζs = rp.ζs # Surface heights

    Nq = length(qs) # Number of wave numbers
    Np = length(ps) # Number of surface points
    Nk = length(ks) # Number of incoming wave numbers
    Nx = length(ζs) # Number of surface points
    Ni = rp.Ni + 1 # Order of surface power expansion

    Npq = zeros(ComplexF64, Np, Nq) # Nₚ⁺(p|q) matrix
    Rqk = Matrix{ComplexF64}(undef, Nq, Nk) # Rₚ(q|k) Solution vector
    Npk = zeros(ComplexF64, Np, Nk) # Nₚ⁻(p|k) vector


    α(q::Float64)::ComplexF64 = √Complex(ε - q^2)

    α0(q::Float64)::ComplexF64 =
        (abs(q) < 1.0) ?
        (√(ε0 - q^2)) :
        (abs(q) > 1.0) ?
        (im * √(q^2 - 1.0)) :
        0.0im

    g_size = Nx ÷ 2 + 1
    gs = Matrix{ComplexF64}(undef, g_size, Ni)
    for n in axes(gs, 2)
        gs[:, n] = rp.FT_plan * ζs .^ (n - 1)
    end



    # Symmetrize gs access -NQ,...,-1, 0, 1,...,NQ
    re_gs(i, n) = (i < 0) ? -gs[-i+1, n] : gs[i+1, n]
    # Fix p-q indexing 
    pq_idx(i, j) = -g_size + i + j - 1

    # Evaluate N⁺ₚ(p,q) matrix
    for n in axes(gs, 2)
        for (i, p) in enumerate(ps)
            for (j, q) in enumerate(qs)
                Npq[i, j] += (p * q + α(p) * α0(q)) * (1.0im)^(n - 1) * (α(p) - α0(q))^(n - 2) / factorial(n - 1) * re_gs(pq_idx(i, j), n)
            end
        end
    end

    # Evaluate N⁻ₚ(p,k) vector
    for n in axes(gs, 2)
        for (i, p) in enumerate(ps)
            for (j, k) in enumerate(ks)
                Npk[i, j] -= (p * k - α(p) * α0(k)) * (1.0im)^(n - 1) * (α(p) + α0(k))^(n - 2) / factorial(n - 1) * re_gs(pq_idx(i, j), n)
            end
        end
    end


    for k in axes(ks, 1)
        Rqk[:, k] = Npq \ Npk[:, k]
    end
    Rqk .*= 2π / rp.Δq

    nans1 = findall(!isfinite, Npq)
    nans2 = findall(!isfinite, Npk)
    nans3 = findall(!isfinite, Rqk)
    if length(nans1) > 0 || length(nans2) > 0 || length(nans3) > 0
        println("Npq:")
        display(Npq[nans1])
        display(nans1)
        println("Npk:")
        display(Npk[nans2])
        display(nans2)
        println("Rqk:")
        display(Rqk[nans3])
        display(nans3)
        throw("Non-finite values in Npq, Npk, or Rqk")
    end
    #= Debug print information
    # println("\ngs:")
    # display(gs)


    println("\nN⁺ₚ(p,q):")
    # display(Npq)

    display(Npq[.!isfinite.(Npq)])
    display(findall(!isfinite, Npq))


    println("\nN⁻ₚ(p,k):")
    # display(Npk)

    display(Npk[.!isfinite.(Npk)])
    display(findall(!isfinite, Npk))

    println("\nRₚ(q|k):")
    display(Rqk)
    =#
    return Rqk
end

rp = RayleighParams();
print(rp)
sol = solve_p(rp);
display(size(sol))

θ(k) = asin(k)

mean_DRC(rp::RayleighParams, sol::Matrix{ComplexF64}) = begin
    # Calculate the mean differential reflection coefficient
    # Factor of ω / c0 removed due to scaling
    retval = similar(sol, Float64)
    for (i, q) in enumerate(rp.qs[rp.qs.>-1.0+1e-3.&&rp.qs.<1.0-1e-3])
        for (j, k) in enumerate(rp.ks)
            retval[i, j] = 1.0 / (rp.Lx) / 2π * cos(θ(q))^2 / cos(θ(k)) * abs2(sol[i, j])
        end
    end
    return retval
end

mean_DRC(rp, sol)


function mean_DRC_incoh(rp::RayleighParams, sol::Matrix{ComplexF64})
    # Calculate the mean differential reflection coefficient for incoherent light
    return 1.0 / (rp.Lx) * rp.ω / (2π * rp.c0)
end

plot(rp.ks, abs2.(sum(sol, dims=1))')

function solve_s(rp::RayleighParams)
    # Solve the reduced Rayleigh equations for s-polarized light
    # TODO: Implement this, currently only p-polarized light is supported
    throw("s-polarized Reduced Rayleigh Method not implemented")

    println("""
    Physical parameters:
    c:  $(rp.c)
    ν:  $(rp.ν) ⟹ κ → μ
    μ₀: $(rp.κ0)
    μ:  $(rp.κ)
    ω:  $(rp.ω)
    θ₀: $(rp.θ₀)

    Sizings:
    Q:  $(rp.Q)
    NQ: $(rp.Nq)
    hQ: $(rp.hq)

    Nx: $(rp.Nx)
    Lx: $(rp.Lx)
    hx: $(rp.hx)

    Ni: $(rp.Ni)
    """)
end