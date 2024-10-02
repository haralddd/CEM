push!(LOAD_PATH, "$(@__DIR__)/../RayleighSolver/")
using RayleighSolver

using ProfileView
using BenchmarkTools


function profile_gaussian_surfacegen()
    rp, sp = default_params_for_surface_testing(GaussianSurface(30.0e-9, 100.0e-9))
    function loop_de_loop(sp, rp)
        for _ in 1:1000
            generate_surface!(sp, rp)
        end
    end
    @btime generate_surface!($sp, $rp)
    ProfileView.@profview loop_de_loop(sp, rp)
end

function profile_rectangular_surfacegen()
    rp, sp = default_params_for_surface_testing(RectangularSurface(30.0e-9, 0.82, 1.97))
    function loop_de_loop(sp, rp)
        for _ in 1:1000
            generate_surface!(sp, rp)
        end
    end
    @btime generate_surface!($sp, $rp)
    ProfileView.@profview loop_de_loop(sp, rp)
end

function profile_solver_components()
    surf = GaussianSurface(30.0e-9, 100.0e-9)
    sp, rp = default_params_for_surface_testing(surf)
    pc = SimPreCompute(rp)
    solve!(sp, rp, pc)
end