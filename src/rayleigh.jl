#= Implements optical scattering under the reduced Rayleigh equations
# Effectively solves Maxwell's equations for a singularly polarized wave
# on a partially symmetric rough surface. Satisfying the boundary conditions

# The reduced Rayleigh equations are a set of coupled integral equations
# which assume that far field scattering conditions (singular direction, up/down in 1D)
# can be used all the way down to the rough surface boundary even though for strongly
# rough surfaces one can get multiple scattering events.
=#

using FFTW

@enum Polarization p s
@enum SurfType flat gauss

function is_power_two(n::Int)::Bool
    # Check if n is a power of 2
    (n & (n - 1)) == 0
end

struct RayleighParams
    ν::Polarization # Polarization [p, s]
    c0 # [m/s], Speed of light in vacuum
    κ0  # ν: ϵ₀ / μ₀
    κ   # Polarization dependent constant
    ω   # Angular frequency
    θ0  # Angle of incidence
    k   # Incidence wave number (k ≡ k∥)

    # Sizings
    Q   # Truncated wave number, q = (-∞,∞) -> q_n = (-Q/2, Q/2)
    Nq  # Number of wave numbers
    Δq  # Wave number spacing [1/m]
    wq  # Weights of the quadrature in q_n

    Nx  # Number of surface points
    Lx  # Surface length [m]
    Δx  # Surface point spacing [m]

    Ni  # Order of surface power expansion

    p   # Scattered wave numbers
    q   # Scattered wave numbers
    ζ   # Surface heights

    FT_plan # Planned Fourier transform of surface points

    RayleighParams(; ν::Polarization=p, surf_t::SurfType=flat, κ0=1.0, κ=1.0, λ=400e-9, θ0=0.0, Q_mult=4, Nq=128, Lx=1.0e-4, Ni=10) = begin
        c0 = 299_792_458.0
        ω = 2π * c0 / λ

        # Assertions and warnings
        @assert Lx / λ > 100.0 "Surface length must be much larger than 1 wavelength, Lx ≫ λ, but is Lx:$Lx and λ:$λ."
        @assert Q_mult > 2 "Q_mult must be greater than 2, but is $Q_mult."
        @assert Nq > 2 "Nq must be greater than 2, but is $Nq."

        # Check if Nq is a power of 2, fits better in cache
        is_power_two(Nq) && @warn("Recommended to use a simple power of 2ⁿ for Nq. Val: $Nq.")

        k = ω / c0 * sin(θ0)

        Q = Q_mult * ω / c0
        Δq = Q / (2 * Nq)

        # Non-inclusive, such that Nq splits the range evenly.
        p = range(-Q / 2 - Δq, Q / 2; length=Nq)
        q = range(-Q / 2, Q / 2 - Δq; length=Nq)

        wq = ones(Nq)
        wq[1] = wq[end] = 3.0 / 8.0
        wq[2] = wq[end-1] = 7.0 / 6.0
        wq[3] = wq[end-2] = 23.0 / 24.0


        if 0 #surf_t == gauss
        # TODO: Tie to Gaussian surface gen. impl.
        else
            ζ = zeros(Nq) # Flat surface
        end

        Nx = 2 * Nq
        Δx = Lx / Nx

        new(ν, c0, κ0, κ, ω, θ0, k, Q, Nq, Δq, wq,
            Nx, Lx, Δx, p, q, ζ, Ni,
            plan_fft(-0.5Lx:Δx:0.5Lx-Δx))
    end
end

function print(rp::RayleighParams)
    names = fieldnames(typeof(rp))
    for name in names
        field = getfield(rp, name)
        println("$name:\t$field")
    end
end

function solve_p(rp::RayleighParams)
    # Solve the reduced Rayleigh equations for p-polarized light

    ε0 = rp.κ0 # Permittivity of free space
    ε = rp.κ # Permittivity of the medium
    ω = rp.ω # Angular frequency


    # Make ranges
    xs = range(-0.5rp.Lx, 0.5rp.Lx - rp.hx; length=rp.Nx) # Coordinates
    ζs = zeros(rp.Nx) # Surface height, flat surface

    qs = range(-0.5Q, 0.5Q - hq; length=rp.Nq) # Wave numbers
    ps = range(-0.5Q, 0.5Q - hq; length=rp.Nq) # Wave numbers

    α(q::ComplexF64)::ComplexF64 = √(ε * (ω / c)^2 - q^2)

    α0(q::ComplexF64)::ComplexF64 =
        (abs(q) < ω / c) ?
        (√(ε0 * (ω / c)^2 - q^2)) :
        (abs(q) > ω / c) ?
        (im * √(q^2 - (ω / c)^2)) :
        0.0im

    # Coordinate space is the set of vectors (p, q), so all operations are in this space,
    #   i.e. f(p) * g(q)' makes the matrix dependent on p and q along each axis


    # For p-polarized light
    Nₚ⁺(p, q) = +(p * q + α(p, ω) * α₀(q, ω)) / (α(p, ω) - α₀(q, ω)) * I(α(p, ω) - α₀(q, ω), p - q)

    Nₚ⁻(p, q) = -(p * q - α(p, ω) * α₀(q, ω)) / (α(p, ω) + α₀(q, ω)) * I(α(p, ω) + α₀(q, ω), p - q)

    # Constrained by u = p - q', I should be a matrix of size N x N




    # Solve the equation ∫dq / 2π N⁺(p, q) R(q, k) = N⁻(p, k), i.e. sum over q

    Npq = Matrix{ComplexF64}(undef, N, N) # Nₚ⁺(p|q) matrix
    Rqk = Matrix{ComplexF64}(undef, N) # Rₚ(q|k) Solution vector
    Npk = Matrix{ComplexF64}(undef, N) # Nₚ⁻(p|k) matrix


    for (n, q) in enumerate(qs)
        for (m, p) in enumerate(ps)
            Npq[m, n] = Nₚ⁺(p, q)
        end
    end

    for (m, p) in enumerate(ps)
        Npk[m, l] = Nₚ⁻(p, k)
    end

    display("Nₚ⁺ = ")
    display(Npq)
    display("Nₚ⁻ = ")
    display(Npk)


    Rqk = Npk / Npq

    display("Rₚ(q|k) = ")
    display(Rqk)


end

solve_p()




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