using RayleighSolver

"For debug purposes, mostly"

function default_config_creation()::SolverData
    return SolverData(Parameters())
end

function default_gaussian_config(iters = 10000, Lx=100*632.8e-9)::SolverData
    return SolverData(Parameters(;Lx=Lx, surf=GaussianSurface(30.0e-9, 100.0e-9)), iters)
end

function default_rectangular_config(iters = 10000)::SolverData
    return SolverData(Parameters(surf=RectangularSurface(30.0e-9, 0.70, 1.30)), iters)
end

function config_glass_isotropic()
    return SolverData(Parameters(below=Isotropic(2.25, 1.0)))
end

function config_silver_isotropic()
    return SolverData(Parameters(below=Isotropic(-7.5 + 0.24im, 1.0)))
end

function config_fresnel_silver(Nθ)
    return SolverData(Parameters(
        θs=range(0.0, 90.0, Nθ),
        Nx=2*2048,
        surf=FlatSurface(),
        below=Isotropic(-7.5 + 0.24im, 1.0)))
end

function config_fresnel_glass(Nθ)
    return SolverData(Parameters(
        θs=range(0.0, 90.0, Nθ),
        Nx=2*2048,
        surf=FlatSurface(),
        below=Isotropic(2.25, 1.0)
    ))
end

function config_default_uniaxial()
    eps_para = 2.0
    eps_perp = 5*eps_para
    return SolverData(Parameters(
        below=Uniaxial(eps_perp, eps_para, 1.0, 1.0)
    ))
end