using RayleighSolver

const eps = -7.5 + 0.24im
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
    below = Uniaxial(eps*ar, eps, 1.0 + 0.0im, 1.0 + 0.0im)
    paramsconf = ParametersConfig(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)
    params = Parameters(paramsconf)

    data = SolverData(params, ens_iters, :reduced)
    solve_ensemble!(data)
    save_solver_data("output/sweep/metallic/$(ar).jld2", data)
end

hyperbolic_ratios = [-0.5, 1.0, -1.0, -1.5]
mkpath("output/sweep/hyperbolic_type1")
for ar in hyperbolic_ratios
    below = Uniaxial(eps, eps*ar, 1.0 + 0.0im, 1.0 + 0.0im)
    paramsconf = ParametersConfig(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)
    params = Parameters(paramsconf)

    data = SolverData(params, ens_iters, :reduced)
    solve_ensemble!(data)
    save_solver_data("output/sweep/hyperbolic_type1/$(ar).jld2", data)
end

mkpath("output/sweep/hyperbolic_type2")
for ar in hyperbolic_ratios
    below = Uniaxial(eps*ar, eps, 1.0 + 0.0im, 1.0 + 0.0im)
    paramsconf = ParametersConfig(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)
    params = Parameters(paramsconf)

    data = SolverData(params, ens_iters, :reduced)
    solve_ensemble!(data)
    save_solver_data("output/sweep/hyperbolic_type2/$(ar).jld2", data)
end
    