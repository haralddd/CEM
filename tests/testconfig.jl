using RayleighSolver

"For debug purposes, mostly"

function default_config_creation()::SolverData
    return SolverData(ParametersConfig())
end

function default_gaussian_config(iters = 10000, Lx=100)::SolverData
    return SolverData(ParametersConfig(;Lx=Lx, surf=GaussianSurface(30.0e-9, 100.0e-9)), iters)
end

function default_rectangular_config(iters = 10000, Lx=100)::SolverData
    return SolverData(ParametersConfig(;Lx=Lx, surf=RectangularSurface(30.0e-9, 0.70, 1.30)), iters)
end

function config_glass_isotropic()
    return SolverData(ParametersConfig(below=Isotropic(2.25, 1.0)))
end

function config_silver_isotropic()
    return SolverData(ParametersConfig(below=Isotropic(-7.5 + 0.24im, 1.0)))
end

function config_fresnel_silver(Nθ)
    return SolverData(ParametersConfig(
        θs=range(0.0, 90.0, Nθ),
        Nx=2*2048,
        surf=FlatSurface(),
        below=Isotropic(-7.5 + 0.24im, 1.0)))
end

function config_fresnel_glass(Nθ)
    return SolverData(ParametersConfig(
        θs=range(0.0, 90.0, Nθ),
        Nx=2*2048,
        surf=FlatSurface(),
        below=Isotropic(2.25, 1.0)
    ))
end

function config_default_uniaxial(;solver_type=:full)
    eps_para = 2.0
    eps_perp = 2.0*eps_para
    return SolverData(ParametersConfig(
        below=Uniaxial(eps_perp, eps_para, 1.0, 1.0)
    ),1,solver_type)
end