push!(LOAD_PATH, "$(@__DIR__)/../RayleighSolver/")
using RayleighSolver

"For debug purposes, mostly"

function default_config_creation()::SolverData
    return SolverData(SimParams())
end

function default_gaussian_config(iters = 10000, Lx=100*632.8e-9)::SolverData
    return SolverData(SimParams(;Lx=Lx, surf=GaussianSurface(30.0e-9, 100.0e-9)), iters)
end

function default_rectangular_config(iters = 10000)::SolverData
    return SolverData(SimParams(surf=RectangularSurface(30.0e-9, 0.70, 1.30)), iters)
end

function config_glass_isotropic()
    return SolverData(SimParams(below=Isotropic(2.25, 1.0)))
end

function config_silver_isotropic()
    return SolverData(SimParams(below=Isotropic(-7.5 + 0.24im, 1.0)))
end

function config_fresnel_silver(Nθ)
    return SolverData(SimParams(
        θs=range(0.0, 90.0, Nθ),
        Nx=2048,
        surf=FlatSurface(),
        below=Isotropic(-7.5 + 0.24im, 1.0)))
end

function config_fresnel_glass(Nθ)
    return SolverData(SimParams(
        θs=range(0.0, 90.0, Nθ),
        Nx=2048,
        surf=FlatSurface(),
        below=Isotropic(2.25, 1.0)
    ))
end