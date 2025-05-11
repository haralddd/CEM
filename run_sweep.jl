using RayleighSolver
using Dates

const eps_re = -7.5
const eps_im = 0.24
const eps = eps_re + eps_im*1.0im
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
sim_name = ARGS[1]
base_path = "output/$sim_name"
mkpath(base_path)
mkpath("$base_path/metallic")
mkpath("$base_path/hyperbolic_type1")
mkpath("$base_path/hyperbolic_type2")

metal_ratios = [0.5, 1.0, 1.5]
for ar in metal_ratios
    eps_perp = eps_re*ar + eps_im*1.0im
    below = Uniaxial(eps_perp, eps, 1.0 + 0.0im, 1.0 + 0.0im)
    paramsconf = ParametersConfig(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)
    params = Parameters(paramsconf)

    data = SolverData(params, ens_iters, :reduced)
    solve_ensemble!(data)
    save_solver_data("$base_path/metallic/$(ar).jld2", data)
end

hyperbolic_ratios = [-0.5, -1.0, -1.5]
for ar in hyperbolic_ratios
    eps_para = eps_re*ar + eps_im*1.0im
    below = Uniaxial(eps, eps_para, 1.0 + 0.0im, 1.0 + 0.0im)
    paramsconf = ParametersConfig(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)
    params = Parameters(paramsconf)

    data = SolverData(params, ens_iters, :reduced)
    solve_ensemble!(data)
    save_solver_data("$base_path/hyperbolic_type1/$(ar).jld2", data)
end

for ar in hyperbolic_ratios
    eps_perp = eps_re*ar + eps_im*1.0im
    below = Uniaxial(eps_perp, eps, 1.0 + 0.0im, 1.0 + 0.0im)
    paramsconf = ParametersConfig(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)
    params = Parameters(paramsconf)

    data = SolverData(params, ens_iters, :reduced)
    solve_ensemble!(data)
    save_solver_data("$base_path/hyperbolic_type2/$(ar).jld2", data)
end
    