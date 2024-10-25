push!(LOAD_PATH, "$(@__DIR__)/../RayleighSolver/")
using RayleighSolver
using Plots
using LinearAlgebra

function test_reciprocity()
    surf = GaussianSurface(30.0e-9, 100.0e-9)
    Q = 4
    Nq = 2048+1
    valid_qs = LinRange(-Q/2, Q/2, Nq)
    valid_ks = valid_qs[-1.0 .< valid_qs .< 1.0]
    @assert all(valid_ks .== .-reverse(valid_ks))
    spa = SimParams(
        lambda=632.8e-9,
        Q=4,
        Nq=Nq,
        ks=valid_ks,
        L=10.0e-6,
        Ni=10,
        surf=surf,
        rescale=true
    )
    sp = SimPrealloc(spa.Nq, length(spa.ks))

    ks = spa.ks
    qs = spa.qs
    Nk = length(ks)
    rev_ks = reverse(ks)
    @assert all(ks .== .-rev_ks)

    @info "generate_surface!"
    @time generate_surface!(sp, spa)

    @info "SimPreCompute"
    @time pc = SimPreCompute(spa)
    validate(pc)

    @time solve!(sp, spa, pc)

    pre(q, k) = √((alpha0(q))/(alpha0(k)))

    Δ = Matrix{Float64}(undef, Nk, Nk)

    rev_idx(idx, N) = N - idx + 1

    for (i, ki) in enumerate(spa.kis)
        for (j, kj) in enumerate(spa.kis)
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
    sp, spa = default_config_creation()
    M = 1
    @time solve_MDRC!(sp, spa, 1)

    M = sp.Mpq
    hm = heatmap(log10.(abs2.(M)), size=(800,800))
    display(hm)
end

function test_solver_surf(surf::T) where T<:RandomSurface
    spa, sp = default_params_for_surface_testing(surf)

    M = 100
    @time qs, coh, incoh = solve_MDRC!(spa, sp, M)

    plt_coh = plot(qs, coh, yscale=:log10, title="$(typeof(surf))\nCoherent MDRC")
    plt_incoh = plot(qs, incoh, yscale=:log10, title="$(typeof(surf))\nIncoherent MDRC")
    return plot(plt_coh, plt_incoh, layout=(1, 2))
end

function test_solver()
    plt1 = test_solver_surf(GaussianSurface(30.0e-9, 100.0e-9))
    plt2 = test_solver_surf(RectangularSurface(30.0e-9, 0.82, 1.97))
    plot(plt1, plt2, layout=(2, 1), size=(800, 800))
end

function test_crystal_precompute()
    ε = 2.25 + 1e-4im
    lambda = 632.8e-9
    Q = 4
    Nq = 2048+1
    ks = [sind(20.0)]
    L = 10.0e-6
    Ni = 3
    surf = GaussianSurface(30.0e-9, 100.0e-9)

    rp_isotropic = SimParams(
        lambda=lambda,
        Q=Q,
        Nq=Nq,
        ks=ks,
        L=L,
        Ni=Ni,
        surf=surf,
        rescale=true,
        above=Vacuum(),
        below=Isotropic(ε, 1.0)
    )
    @time pc_i = SimPreCompute(rp_isotropic)
    validate(pc_i)

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
    @time pc_c = SimPreCompute(rp_crystal)
    validate(pc_c)


    ln_Ni = plot()
    ln_Nc = plot()
    ln_Ndiff = plot()

    hm_Mi = surface(title="Isotropic")
    hm_Mc = surface(title="Uniaxial Crystal")
    hm_Mdiff = surface(title="Error")

    Ni = abs2.(pc_i.Npkn)
    Nc = abs2.(pc_c.Npkn)

    Mi = log10.(abs2.(pc_i.Mpqn))
    Mc = log10.(abs2.(pc_c.Mpqn))

    Ndiff = abs2.(pc_i.Npkn .- pc_c.Npkn)
    Mdiff = abs2.(pc_i.Mpqn .- pc_c.Mpqn)

    for n in axes(pc_i.Mpqn, 3)
        plot!(ln_Ni, Ni[:,1,n], label = "Isotropic N, n=$n")
        plot!(ln_Nc, Nc[:,1,n], label = "Uniaxial Crystal N, n=$n")
        plot!(ln_Ndiff, Ndiff[:,1,n], label = "Error N, n=$n")

        surface!(hm_Mi, Mi[:,:,n] .+ n)
        surface!(hm_Mc, Mc[:,:,n] .+ n)
        surface!(hm_Mdiff, Mdiff[:,:,n] .+ n)
    end
    

    plot(
        hm_Mi, ln_Ni,
        hm_Mc, ln_Nc,
        hm_Mdiff, ln_Ndiff,
        layout=(3, 2), size=(800, 800)
        ) |> display
    
    @show maximum(abs2.(pc_i.Mpqn) .- abs2.(pc_c.Mpqn))
    @show maximum(Mdiff)
    @show maximum(Mdiff[:,:,1])
    @show maximum(Mdiff[:,:,2])
    @show maximum(Mdiff[:,:,3])
    @show argmax(Mdiff)
    @show maximum(abs2.(pc_i.Npkn) .- abs2.(pc_c.Npkn))
    @show maximum(Ndiff)
    @show argmax(Ndiff)
    nothing
end
