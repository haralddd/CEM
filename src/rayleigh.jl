#= Implements optical scattering under the reduced Rayleigh equations
# Effectively solves Maxwell's equations for a singularly polarized wave
# on a partially symmetric rough surface. Satisfying the boundary conditions

# The reduced Rayleigh equations are a set of coupled integral equations
# which assume that far field scattering conditions (singular direction, up/down in 1D)
# can be used all the way down to the rough surface boundary even though for strongly
# rough surfaces one can get multiple scattering events.
=#

using FFTW



😄(🐱) = '🐶'

😄('🐴')


🐴 = 🔥 -> 🔥^2 * 😄(🔥)^4

🐴("🐈")

for j in 1:10
    str = ""
    for i in 1:10
        str *= Char(0x1F430 + i + j * 10)
    end
    println(str)
end


L = 1.0 # Length of surface
N = 1000 # Number of surface points
Δξ = L / N # Surface point spacing

x = -0.5L:Δξ:0.5L-Δξ # Surface points

Fx = plan_fft(x) # Plan Fourier transform of surface points 
ζ = zeros(N) # Surface height, flat surface

Iq(γ, n) = Fx * exp(-im * γ * ζ[n]) # I(γ | q) = ∫dx exp(-iγζ(x)) ⋅ exp(-iqx)

M⁺(p, q) = +(
    (p + κ(ω) * q) * (p - q) / (α(p, ω))
)