push!(LOAD_PATH, "$(@__DIR__)/../RayleighSolver/")
using RayleighSolver
using Plots

function test_solve(surf::T) where T<:SurfaceParams
    rp, sp = default_params_for_surface_testing(surf)

    M = 10
    qs, coh, incoh = solve_MDRC!(rp, sp, M)
    plt_coh = plot(qs, coh, yscale=:log10, title="$surf coherent MDRC")
    plt_incoh = plot(qs, incoh, yscale=:log10, title="$surf incoherent MDRC")
    return plot(plt_coh, plt_incoh, layout=(1, 2))
end



function test_solver()
    plt1 = test_solve(GaussianSurfaceParams(30.0e-9, 100.0e-9))
    plt2 = test_solve(RectSurfaceParams(30.0e-9, 0.82, 1.97))
    plot(plt1, plt2, layout=(2, 1), size=(800, 800))
end
