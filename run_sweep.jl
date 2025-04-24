using RayleighSolver
using CairoMakie
using LaTeXStrings
using JLD2, FileIO

const eps_perp = -7.5 + 0.24im
const surf_rms = 10.0e-9
const km = 0.782
const kp = 1.366
const θs = [10.0]
const λ = 457.9e-9
const Nx = 4096
const Lx = 200
const Ni = 10
const ens_iters = 10000

surface = RectangularSurface(surf_rms, km, kp)
above = Vacuum()

# Parameter ranges to sweep

metal_ratios = [0.5, 1.0, 1.5]
mkpath("output/sweep/metallic")
for ar in metal_ratios
    below = Uniaxial(eps_perp, eps_perp*ar, 1.0 + 0.0im, 1.0 + 0.0im)
    params = Parameters(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)

    data = SolverData(params, ens_iters, :reduced)
    solve_ensemble!(data)
    save("output/sweep/metallic/$(ar).jld2", data)
end

hyperbolic_ratios = [-0.5, 1.0, -1.0, -1.5]
mkpath("output/sweep/hyperbolic")
for ar in hyperbolic_ratios
    below = Uniaxial(eps_perp, eps_perp*ar, 1.0 + 0.0im, 1.0 + 0.0im)
    params = Parameters(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)

    data = SolverData(params, ens_iters, :reduced)
    solve_ensemble!(data)
    save("output/sweep/hyperbolic/$(ar).jld2", data)
end