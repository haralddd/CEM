push!(LOAD_PATH, "$(@__DIR__)/../RayleighSolver/")
using RayleighSolver
using Plots
using LinearAlgebra


function test_reciprocity()
    surf = GaussianSurfaceParams(30.0e-9, 100.0e-9)
    rp = RayleighParams(
        nu=p,
        eps=ComplexF64(2.25),
        mu=ComplexF64(1.0),
        lambda=632.8e-9,
        Q=4,
        Nq=1024,
        ks=[sind(-20.0), sind(-10.0), sind(10.0), sind(20.0)],
        L=10.0e-6,
        Ni=10,
        surf=surf,
        rescale=true
    )
    sp = SimulationPreAlloc(rp.Nq, length(rp.ks))
    generate!(sp, rp)
    @time Mpk_pre, Npk_pre = precalc(rp)

    @assert all(isfinite.(Mpk_pre))
    @assert all(isfinite.(Npk_pre))

    solve!(sp, rp, Mpk_pre, Npk_pre)

    pre(q, k) = √((alpha0(q))/(alpha0(k)))
    a1 = pre.(rp.qs, rp.ks[1]) .* sp.Npk[:, 1]
    a2 = pre.(rp.qs, rp.ks[2]) .* sp.Npk[:, 2]
    a3 = pre.(rp.qs, rp.ks[3]) .* sp.Npk[:, 3]
    a4 = pre.(rp.qs, rp.ks[4]) .* sp.Npk[:, 4]

    plt1 = plot(log10.(abs.(real.(a1) .- real.(a4))), label="log real Δ, θ = -20, 20")
    plot!(log10.(abs.(imag.(a1) .- imag.(a4))), label="log imag Δ, θ = -20, 20")

    plt2 = plot(log10.(abs.(real.(a2) .- real.(a3))), label="log real Δ, θ = -10, 10")
    plot!(log10.(abs.(imag.(a2) .- imag.(a3))), label="log imag Δ, θ = -10, 10")


    plot(plt1, plt2, layout=(2, 1), size=(800, 800))
end

function test_solve(surf::T) where T<:SurfaceParams
    rp, sp = default_params_for_surface_testing(surf)

    M = 10
    qs, coh, incoh = solve_MDRC!(rp, sp, M)

    # Should satisfy reciprocity


    hm = heatmap(log10.(abs.(real.(sp.Mpq))))
    hm2 = heatmap(log10.(abs.(imag.(sp.Mpq))))
    plt_coh = plot(qs, coh, yscale=:log10, title="$surf coherent MDRC")
    plt_incoh = plot(qs, incoh, yscale=:log10, title="$surf incoherent MDRC")
    return plot(plt_coh, plt_incoh, hm, hm2, layout=(2, 2))
end



function test_solver()
    plt1 = test_solve(GaussianSurfaceParams(30.0e-9, 100.0e-9))
    plot(plt1, size=(800, 800))

    # plt2 = test_solve(RectSurfaceParams(30.0e-9, 0.82, 1.97))
    # plot(plt1, plt2, layout=(2, 1), size=(800, 800))
end
