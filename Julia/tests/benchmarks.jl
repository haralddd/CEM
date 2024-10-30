push!(LOAD_PATH, "$(@__DIR__)/../RayleighSolver/")
using RayleighSolver

using ProfileView
using BenchmarkTools


function profile_gaussian_surfacegen()
    spa, sp = default_params_for_surface_testing(GaussianSurface(30.0e-9, 100.0e-9))
    function loop_de_loop(sp, spa)
        for _ in 1:1000
            generate_surface!(sp, spa)
        end
    end
    @btime generate_surface!($sp, $spa)
    ProfileView.@profview loop_de_loop(sp, spa)
end

function profile_rectangular_surfacegen()
    spa, sp = default_params_for_surface_testing(RectangularSurface(30.0e-9, 0.82, 1.97))
    function loop_de_loop(sp, spa)
        for _ in 1:1000
            generate_surface!(sp, spa)
        end
    end
    @btime generate_surface!($sp, $spa)
    ProfileView.@profview loop_de_loop(sp, spa)
end

function profile_solver_components()
    surf = GaussianSurface(30.0e-9, 100.0e-9)
    sp, spa = default_params_for_surface_testing(surf)
    pc = SimPreCompute(spa)
    solve_single!(sp, spa, pc)
end

function profile_crystal_precompute()
    ε = 2.25 + 1e-4im
    lambda = 632.8e-9
    Q = 4
    Nq = 512 + 1
    ks = [sind(20.0)]
    L = 10.0e-6
    Ni = 3
    surf = GaussianSurface(30.0e-9, 100.0e-9)
    rp_crystal = SimParams(
        lambda=lambda,
        Q=Q,
        Nq=Nq,
        ks=ks,
        L=L,
        Ni=Ni,
        surf=surf,
        rescale=true,
        above=UniaxialCrystal(1.0, 1.0, 1.0, 1.0),
        below=UniaxialCrystal(ε, ε, 1.0, 1.0)
    )

    @benchmark SimPreCompute($rp_crystal)
end