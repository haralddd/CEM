push!(LOAD_PATH, "$(@__DIR__)/../RayleighSolver/")
using RayleighSolver
using Plots
using LinearAlgebra


function test_reciprocity()
    surf = GaussianSurfaceParams(30.0e-9, 100.0e-9)
    Q = 4
    Nq = 2*4096+1
    valid_qs = LinRange(-Q/2, Q/2, Nq)
    mid = Nq÷2 + 1
    ks = [valid_qs[mid - mid÷3], valid_qs[mid - mid÷4], valid_qs[mid + mid÷4], valid_qs[mid + mid÷3]]
    display(ks)
    @assert ks[1] == -ks[end]
    @assert ks[2] == -ks[end-1]
    rp = RayleighParams(
        nu=p,
        eps=ComplexF64(2.25),
        mu=ComplexF64(1.0),
        lambda=632.8e-9,
        Q=4,
        Nq=Nq,
        ks=ks,
        L=10.0e-6,
        Ni=10,
        surf=surf,
        rescale=true
    )
    display(rp.ks)
    sp = SimulationPreAlloc(rp.Nq, length(rp.ks))
    generate!(sp, rp)
    @time Mpk_pre, Npk_pre = precalc(rp)

    @assert all(isfinite.(Mpk_pre))
    @assert all(isfinite.(Npk_pre))

    solve!(sp, rp, Mpk_pre, Npk_pre)

    pre(q, k) = √((alpha0(q))/(alpha0(k)))
    a1 = pre.(reverse(rp.qs), rp.ks[1]) .* reverse(sp.Npk[:, 1])
    a2 = pre.(reverse(rp.qs), rp.ks[2]) .* reverse(sp.Npk[:, 2])
    a3 = pre.(rp.qs, rp.ks[3]) .* sp.Npk[:, 3]
    a4 = pre.(rp.qs, rp.ks[4]) .* sp.Npk[:, 4]


    plt1 = plot(log10.(abs.(a1 .- a4)), label="error log10, k = $(ks[1]), $(ks[4])")

    plt2 = plot(log10.(abs.(a2 .- a3)), label="error log10, S($(ks[1])|-q) - S($(ks[4])|q)")


    plt = plot(plt1, plt2, layout=(2, 1), size=(800, 800))
    display(plt)
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
