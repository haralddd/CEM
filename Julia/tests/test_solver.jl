include("testconfig.jl")
# using Plots
using LinearAlgebra
using Statistics

function test_observation()
    N = 10000000
    A = rand(N)*1.e-6 .+ 1.0

    @info "mean(A) = $(mean(A))"

    observable = 0.0
    for n in eachindex(A)
        observable = observe(observable, A[n], n)
    end

    @info "observable mean: A = $observable"
    @info "difference: $(observable - mean(A))"
end

function test_reciprocity()
    surf = GaussianSurface(30.0e-9, 100.0e-9)
    Q = 4
    Nq = 2048

    dq = Q / (Nq - 1)
    qs = -Q/2:dq:Q/2
    ks = qs[-1.0 .< qs .< 1.0]
    Nk = length(ks)

    spa = SimParams(
        lambda=632.8e-9,
        Q=Q,
        Nq=Nq,
        ks=ks,
        Lx=10.0e-6,
        Ni=10,
        surf=surf,
        rescale=true
    )
    data = SolverData(spa, 1)

    @info "SimPreCompute"
    @time precompute!(data)
    spa = data.spa
    sp = data.sp
    pd = sp.p_data
    sd = sp.s_data


    @info "generate_surface!"
    @time generate_surface!(sp, spa)

    @info "solve_single!"
    @time solve_single!(data)

    Δ_p = fill(1e-20, Nk, Nk)
    Δ_s = fill(1e-20, Nk, Nk)

    rev_idx(idx) = Nk - idx + 1

    for (j, kj) in enumerate(spa.kis)
        for (i, qi) in enumerate(spa.kis)

            k = qs[kj]
            q = qs[qi]

            mi = rev_idx(i)
            mqj = spa.rev_kis[j]
            mqi = spa.rev_kis[i]

            mk = qs[mqj]
            mq = qs[mqi]

            @assert ks[mi] == -ks[i]
            @assert mk == -k
            @assert mq == -q

            a1 = alpha0(q)
            a2 = alpha0(k)

            b1 = alpha0(mk)
            b2 = alpha0(mq)

            if a2 ≈ 0.0 || b2 ≈ 0.0 continue end

            Sqk_p = sqrt(a1/a2) * pd.Npk[qi, j]
            Skq_p = sqrt(b1/b2) * pd.Npk[mqj, mi]

            Sqk_s = sqrt(a1 / a2) * sd.Npk[qi, j]
            Skq_s = sqrt(b1 / b2) * sd.Npk[mqj, mi]

            Δ_p[i, j] = abs(Sqk_p - Skq_p)
            Δ_s[i, j] = abs(Sqk_s - Skq_s)
        end
    end

    hm_p = heatmap(ks, ks, log10.(Δ_p), size=(800, 800))
    hm_s = heatmap(ks, ks, log10.(Δ_s), size=(800, 800))
    display(plot(hm_p, hm_s))
    display("Reciprocity in P polarization, maximum error: $(maximum(Δ_p))")
    display("Reciprocity in S polarization, maximum error: $(maximum(Δ_s))")
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
