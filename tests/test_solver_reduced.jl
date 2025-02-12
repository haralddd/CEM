include("testconfig.jl")
using Plots
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


function plot_singlesolve()
    data = default_config_creation()
    @info "precompute"
    @time precompute!(data)

    alloc = Preallocated(data.params)

    @info "surface generator"
    @time generate_surface!(alloc, data.params)

    @info "solve single"
    @time solve_single!(alloc, data)

    qs = data.params.qs
    ks = data.params.ks
    mask = qs .>= -1.0 .&& qs .<= 1.0
    Npk = alloc.PNpk
    # heatmap(log10.(abs2.(data.sp.p_data.Mpqn[:,:,7])))
    # display(data.params.qs)
    plot()
    plot!(qs[mask], [log10.(abs2.(Npk[mask, i])) for i in axes(Npk,2)])
    vline!(ks)
    # heatmap(data.params.qs, data.params.ks, )
end

plot_singlesolve()

function plot_Rmean()
    data = SolverData(Parameters(), 100)
    solve_MDRC!(data)

    qs = data.params.qs
    ks = data.params.ks
    mask = qs .>= -1.0 .&& qs .<= 1.0
    R = data.P_res.R

    plot()
    plot!(qs[mask], [abs2.(R[mask, i]) for i in axes(R, 2)])
    vline!(ks)
end
plot_Rmean()

function test_save_data()

    data = SolverData(Parameters(), 100)
    solve_MDRC!(data)

    qs = data.params.qs
    ks = data.params.ks
    mask = qs .>= -1.0 .&& qs .<= 1.0
    display(mask)
    display(size(mask))
    R = data.P_res.R

    plt1 = plot()
    plot!(qs[mask], [abs2.(R[mask, i]) for i in axes(R, 2)])
    vline!(ks)

    # rm("test_save.jdl2")
    save_solver_data("test_save", data)
    loaded_data = load_solver_data("test_save")
    R_ld = loaded_data.P_res.R

    plt2 = plot()
    plot!(qs[mask], [abs2.(R_ld[mask, i]) for i in axes(R, 2)])
    vline!(ks)

    plot(plt1, plt2) |> display
    rm("test_save.jld2")
end

test_save_data()


function test_crystal_precompute()
    ε = 2.25 + 1e-4im
    lambda = 632.8e-9
    Q = 4
    Nq = 2048+1
    ks = [sind(20.0)]
    L = 10.0e-6
    Ni = 3
    surf = GaussianSurface(30.0e-9, 100.0e-9)

    rp_isotropic = Parameters(
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

    rp_crystal = Parameters(
        lambda=lambda,
        Q=Q,
        Nq=Nq,
        ks=ks,
        L=L,
        Ni=Ni,
        surf=surf,
        rescale=true,
        above=Uniaxial(1.0, 1.0, 1.0, 1.0),
        below=Uniaxial(ε, ε, 1.0, 1.0)
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
