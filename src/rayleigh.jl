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
    L = 1.0 # Length of surface
    N = 10 # Number of surface points
    Δξ = L / N # Surface point spacing

    ξs = -0.5L:Δξ:0.5L-Δξ # Surface points

    ζ = zeros(N) # Surface height, flat surface
    qs = fftfreq(N) |> fftshift
    ps = fftfreq(N) |> fftshift
    ks = fftfreq(N) |> fftshift
    ω = 1.0 # Angular frequency

    PF = plan_fft(ξs) # Plan Fourier transform of surface points 



    ε(ω) = 1.0 # Permittivity
    μ(ω) = 0.5 # Permeability
    c = 1.0 # Speed of light

    α(q, ω) = √(ε(ω) * μ(ω) * (ω / c)^2 - q^2)

    α₀(q, ω) =
        (abs(q) < ω / c) ?
        (√(ε(ω) * μ(ω) * (ω / c)^2 - q^2)) :
        (im * √(q^2 - (ω / c)^2)) # In vacuum, ε = μ = 1

    Iq(γ) = PF * exp.(-im * γ .* ζ) # I(γ | q) = ∫dx exp(-iγζ(x)) ⋅ exp(-iqx)

    # Coordinate space is the set of vectors (p, q), so all operations are in this space,
    #   i.e. f(p) * g(q)' makes the matrix dependent on p and q along each axis

    N⁺(p, q) = +(p * q' + α(p, ω) * α₀(q, ω)') / (α(p, ω) - α₀(q, ω)') * Iq(α(p, ω) - α₀(q, ω), p - qx)

    N⁻(p, q) = -(p * q' - α(p, ω) * α₀(q, ω)') / (α(p, ω) + α₀(q, ω)') * Iq(α(p, ω) + α₀(q, ω), p - q)


    Rp = Vector{ComplexF64}(undef, N) # Solution vector, varying in q
    Npq = Matrix{ComplexF64}(undef, N, N) # N⁺ matrix, varying in p and q
    Npk = Matrix{ComplexF64}(undef, N, N) # N⁻ matrix, varying in p and k

    # Solve the equation ∫dq / 2π N⁺(p, q) R(q, k) = N⁻(p, k), i.e. sum over q
    for (n, q) in enumerate(qs)
        for (m, p) in enumerate(ps)
            for (l, k) in enumerate(ks)
                Rp[l] = Npk[p, k] * R[k]
            end
            Npq[m, n] = N⁺(p, q)
            Npk[p, q] = N⁻(p, q)
        end
    end


end

solve()