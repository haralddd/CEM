push!(LOAD_PATH, "$(@__DIR__)/../RayleighSolver/")
using RayleighSolver
using Plots
using LinearAlgebra
using ProfileView
using BenchmarkTools

function test_reciprocity()
    surf = GaussianSurface(30.0e-9, 100.0e-9)
    Q = 4
    Nq = 2048+1
    valid_qs = LinRange(-Q/2, Q/2, Nq)
    valid_ks = valid_qs[-1.0 .< valid_qs .< 1.0]
    @assert all(valid_ks .== .-reverse(valid_ks))
    rp = SimParams(
        lambda=632.8e-9,
        Q=4,
        Nq=Nq,
        ks=valid_ks,
        L=10.0e-6,
        Ni=10,
        surf=surf,
        rescale=true
    )
    sp = SimPrealloc(rp.Nq, length(rp.ks))

    ks = rp.ks
    qs = rp.qs
    Nk = length(ks)
    rev_ks = reverse(ks)
    @assert all(ks .== .-rev_ks)

    @info "generate_surface!"
    @time generate_surface!(sp, rp)

    @info "SimPreCompute"
    @time pc = SimPreCompute(rp)
    validate(pc)

    @time solve!(sp, rp, pc)

    pre(q, k) = √((alpha0(q))/(alpha0(k)))

    Δ = Matrix{Float64}(undef, Nk, Nk)

    rev_idx(idx, N) = N - idx + 1

    for (i, ki) in enumerate(rp.kis)
        for (j, kj) in enumerate(rp.kis)
            q = ks[i]
            k = ks[j]
            i_rev = rev_idx(i, Nk)
            kj_rev = rev_idx(kj, Nq)
            @assert q == -ks[i_rev] "q=$(q) != -ks[i_rev]=$(-ks[i_rev])"
            @assert k == -qs[kj_rev] "k=$(k) != -qs[kj_rev]=$(-qs[kj_rev])"

            Sqk = pre(q, k) * sp.Npk[ki, j]
            Skq = pre(-k, -q) * sp.Npk[kj_rev, i_rev]

            Δ[i, j] = abs(Sqk - Skq)
        end
    end

    hm = heatmap(ks, ks, log10.(Δ), size=(800, 800))
    display(hm)
    display("Reciprocity, maximum error: $(maximum(Δ))")
end

function test_symmetry_isotropic()
    sp, rp = default_config_creation()
    M = 1
    @time solve_MDRC!(sp, rp, 1)

    M = sp.Mpq
    hm = heatmap(log10.(abs2.(M)), size=(800,800))
    display(hm)
end

function test_solver_surf(surf::T) where T<:RandomSurface
    rp, sp = default_params_for_surface_testing(surf)

    M = 100
    @time qs, coh, incoh = solve_MDRC!(rp, sp, M)

    plt_coh = plot(qs, coh, yscale=:log10, title="$(typeof(surf))\nCoherent MDRC")
    plt_incoh = plot(qs, incoh, yscale=:log10, title="$(typeof(surf))\nIncoherent MDRC")
    return plot(plt_coh, plt_incoh, layout=(1, 2))
end

function test_solver()
    plt1 = test_solver_surf(GaussianSurface(30.0e-9, 100.0e-9))
    plt2 = test_solver_surf(RectangularSurface(30.0e-9, 0.82, 1.97))
    plot(plt1, plt2, layout=(2, 1), size=(800, 800))
end

function profile_solver_components()
    surf = GaussianSurface(30.0e-9, 100.0e-9)
    sp, rp = default_params_for_surface_testing(surf)
    pc = SimPreCompute(rp)
    solve!(sp, rp, pc)
end
