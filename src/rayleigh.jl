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



function solve()

    ε₀ = 1.0 # Relative permittivity of vacuum
    μ₀ = 1.0 # Relative permeability of vacuum

    μ = 2.0 # Permeability
    ε = 1.0 # Permittivity



    # Speed of light in vacuum
    c = 299_792_458.0 # [m/s]

    α(q::ComplexF64, ω::Float64)::ComplexF64 = √(ε * μ * (ω / c)^2 - q^2)

    α₀(q::ComplexF64, ω::Float64)::ComplexF64 =
        (abs(q) < ω / c) ?
        (√(ε₀ * μ₀ * (ω / c)^2 - q^2)) :
        (abs(q) > ω / c) ?
        (im * √(q^2 - (ω / c)^2)) :
        0.0im # In vacuum, ε = μ = 1

    L = 1.0e-6 # Length of surface [m]
    Q = 6.0e-6 # Truncated wave number [1/m]
    Nq = 128 # Wave number grid points
    Nx = Nq^2 # Surface grid points, must be able to resolve p-q in Fourier space
    Δξ = L / N # Surface point spacing

    θ₀ = asin(√() / √(ε₀ * μ₀)) # Angle of incidence
    k = ω / c * sin(θ₀)

    ξs = range(-0.5L, 0.5L, Nq) # Surface points
    ζs = zeros(Nx) .|> ComplexF64 # Surface height, flat surface
    qs = fftfreq(N, 2π * Δξ) |> fftshift .|> ComplexF64
    ps = fftfreq(N, 2π * Δξ) .+ 2.0e-8 |> fftshift .|> ComplexF64


    ω = 1.0e6 # Angular frequency

    PF = plan_fft(ζs) # Plan Fourier transform of surface points

    # I(γ|q) = ∫dx exp(-iγζ(x)) ⋅ exp(-iux)
    # Where u is the resulting coordinate space
    I(γ, q) = reduce(+, exp.(-im * γ * ζs) .* exp.(-im * q * ξs))

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