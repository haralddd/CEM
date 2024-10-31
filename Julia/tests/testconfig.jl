push!(LOAD_PATH, "$(@__DIR__)/../RayleighSolver/")
using RayleighSolver

"For debug purposes, mostly"
function config_glass_isotropic()
    ε = 2.25
    lambda = 632.8e-9
    Q = 4
    Nq = 1024
    ks = sind.(0.0:10.0:20.0)
    Lx = 10.0e-6
    Ni = 10
    surf = GaussianSurface(10.0e-9, 50.0e-9)
    spa = SimParams(
        lambda=lambda,
        Q=Q,
        Nq=Nq,
        ks=ks,
        Lx=Lx,
        Ni=Ni,
        surf=surf,
        rescale=true,
        above=Vacuum(),
        below=Isotropic(ε, 1.0)
    )
    return spa
end

function config_silver_isotropic()
    ε = -7.5 + 0.24im
    lambda = 632.8e-9
    Q = 4
    Nq = 1024
    ks = sind.(0.0:10.0:20.0)
    Lx = 10.0e-6
    Ni = 10
    surf = GaussianSurface(10.0e-9, 50.0e-9)
    spa = SimParams(
        lambda=lambda,
        Q=Q,
        Nq=Nq,
        ks=ks,
        Lx=Lx,
        Ni=Ni,
        surf=surf,
        rescale=true,
        above=Vacuum(),
        below=Isotropic(ε, 1.0)
    )
    return spa
end

function default_config_creation()::SolverData
    spa = SimParams(
        lambda=632.8e-9,
        Q=4,
        Nq=1024,
        ks=[sind(10.0), sind(20.0), sind(30.0)],
        Lx=10.0e-6,
        Ni=10,
        surf=GaussianSurface(30.0e-9, 100.0e-9),
        rescale=true
    )

    return SolverData(spa)
end

function default_params_for_surface_testing(surf::T)::Tuple{SimPrealloc,SimParams} where {T<:RandomSurface}
    spa = SimParams(
        lambda=632.8e-9,
        Q=4,
        Nq=1024,
        ks=[sind(10.0), sind(20.0), sind(30.0)],
        Lx=10.0e-6,
        Ni=10,
        surf=surf,
        rescale=true
    )

    sp = SimPrealloc(spa.Nq, length(spa.ks))

    return sp, spa

end