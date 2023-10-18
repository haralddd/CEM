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
using SparseArrays
using BenchmarkTools

struct RayleighParams
    ν = 'p' # Polarization [p, s]

    # Input parameteres
    c0 = 299_792_458.0 # [m/s], Speed of light in vacuum

    κ0 = 1.0 # ν: ϵ₀ / μ₀
    κ = 0.1 # Polarization dependent constant
    ω = 1.0e6 # Angular frequency
    θ₀ = 0.0 # Angle of incidence
    k = ω / c * sin(θ₀)

    # Sizings
    Q = 2 * ω / c # Truncated wave number [1/m]
    NQ = 128 # Number of wave numbers
    ΔQ = Q / NQ # Wave number spacing [1/m]

    Nx = 2 * Nq
    Lx = 1.0e-6 # Surface length [m]
    Δx = Lx / Nx # Surface point spacing [m]

    Ni = 20 # Order of surface expansion
end

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

function solve_p(rp::RayleighParams)
    # Solve the reduced Rayleigh equations for p-polarized light

    println("""
    Physical parameters:
        c:  $(rp.c)
        ν:  $(rp.ν) ⟹ κ → ε
        ε₀: $(rp.κ0)
        ε:  $(rp.κ)
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



    FT = plan_fft(xs)   # Plan Fourier transform of surface points,
    # corresponding to the total space spanned by p-q

    gs = Vector{Matrix{ComplexF64}}(undef, rp.Ni) # Surface expansion Fourier transforms
    for i in eachindex(gs)
        gs[i] = FT * ζs .^ (i - 1)
    end

    gsi = Vector{Matrix{ComplexF64}}(undef, rp.Ni) # Inverse surface expansion Fourier transforms
    for i in eachindex(gsi)
        gsi[i] = FT \ ζs .^ (i - 1)
    end

    # Coordinate space is the set of vectors (p, q), so all operations are in this space,
    #   i.e. f(p) * g(q)' makes the matrix dependent on p and q along each axis


    # For p-polarized light
    Nₚ⁺(p, q) = +(p * q + α(p, ω) * α₀(q, ω)) / (α(p, ω) - α₀(q, ω)) * I(α(p, ω) - α₀(q, ω), p - q)

    Nₚ⁻(p, q) = -(p * q - α(p, ω) * α₀(q, ω)) / (α(p, ω) + α₀(q, ω)) * I(α(p, ω) + α₀(q, ω), p - q)

    # Constrained by u = p - q', I should be a matrix of size N x N




    # Solve the equation ∫dq / 2π N⁺(p, q) R(q, k) = N⁻(p, k), i.e. sum over q

    Npq = Matrix{ComplexF64}(undef, N, N) # Nₚ⁺(p|q) matrix
    Rqk = Matrix{ComplexF64}(undef, N) # Rₚ(q|k) Solution matrix
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

solve()