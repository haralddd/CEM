include("testconfig.jl")
using CairoMakie

function test_reciprocity()
    λ = 632.8e-9
    Nx = 2048
    _Lx = 100*λ
    
    Lx = _Lx * 2π / λ
    dx = Lx / Nx
    dq = 2π/Lx
    Q = π/dx

    qs = -Q/2:dq:Q/2
    mask = -1.0 .< qs .< 1.0
    ks = qs[mask]
    Nk = length(ks)
    θs = asind.(ks)

    @info "λ: $λ"
    @info "Lx: $Lx"
    @info "Nx: $Nx"
    @info "ks: $(ks[1]):$(ks[end])"
    @info "qs has zero: $(any(qs .== 0.0))"


    # below = Isotropic(2.0, 1.0)
    below = Isotropic(-10.0, 1.0)
    spa = SimParams(
        lambda=λ,
        Nx=Nx,
        θs=θs,
        Lx=_Lx,
        below=below
    )
    data = SolverData(spa, 1)

    @info "SimPreCompute"
    @time precompute!(data)


    @info "generate_surface!"
    @time generate_surface!(data.sp, data.spa)

    @info "solve_single!"
    @time solve_single!(data)
    pd = data.sp.p_data
    sd = data.sp.s_data

    _l = lines(θs, abs2.(pd.Npk[mask, 1]), label="θ = $(θs[1])")
    lines!(θs, abs2.(pd.Npk[mask, Nk÷2+1]), label="θ = $(θs[Nk÷2+1])")
    lines!(θs, abs2.(pd.Npk[mask, end]), label="θ = $(θs[end])")
    axislegend()
    display(_l)

    Δ_p = fill(1e-20, Nk, Nk)
    Δ_s = fill(1e-20, Nk, Nk)

    rev_idx(idx) = Nk - idx + 1
    kis = spa.kis
    rev_qis = reverse(eachindex(qs))
    rev_kis = rev_qis[kis]

    @assert(any(pd.Npk .!= 0.0))
    @assert(any(sd.Npk .!= 0.0))

    for (j, kj) in enumerate(kis)
        for (i, qi) in enumerate(kis)

            k = qs[kj]
            q = qs[qi]

            mi = rev_idx(i)
            mqj = rev_kis[j]
            mqi = rev_kis[i]

            mk = qs[mqj]
            mq = qs[mqi]

            @assert ks[mi] == -ks[i]
            @assert mk == -k
            @assert mq == -q

            a1 = alpha0(q)
            a2 = alpha0(k)

            b1 = alpha0(mk)
            b2 = alpha0(mq)

            # if a2 ≈ 0.0 || b2 ≈ 0.0 continue end
            # if a1 ≈ 0.0 || b1 ≈ 0.0 @warn "zero" end

            Sqk_p = sqrt(a1/a2) * pd.Npk[qi, j]
            Skq_p = sqrt(b1/b2) * pd.Npk[mqj, mi]

            Sqk_s = sqrt(a1 / a2) * sd.Npk[qi, j]
            Skq_s = sqrt(b1 / b2) * sd.Npk[mqj, mi]

            Δ_p[i, j] = abs(Sqk_p - Skq_p)
            Δ_s[i, j] = abs(Sqk_s - Skq_s)
        end
    end

    fig = Figure(; size=(800, 400))
    ax, hm1 = heatmap(fig[1, 1], ks, ks, Δ_p)
    Colorbar(fig[1, 2], hm1)
    ax, hm2 = heatmap(fig[1, 3], ks, ks, Δ_s)
    Colorbar(fig[1, 4], hm2)
    display(fig)

    maxidx_p = argmax(Δ_p)[2]
    maxidx_s = argmax(Δ_s)[2]
    qmaxidx_p = kis[maxidx_p]
    qmaxidx_s = kis[maxidx_s]
    rel_err_p = Δ_p[maxidx_p] / abs(pd.Npk[qmaxidx_p, maxidx_p])
    rel_err_s = Δ_s[maxidx_s] / abs(sd.Npk[qmaxidx_s, maxidx_s])

    @info "Reciprocity in P polarization, maximum abs error: $(maximum(Δ_p))"
    @info "Reciprocity in S polarization, maximum abs error: $(maximum(Δ_s))"
    @info "Reciprocity in P polarization, maximum rel error: $(rel_err_p)"
    @info "Reciprocity in S polarization, maximum rel error: $(rel_err_s)"
end

if (abspath(PROGRAM_FILE) == @__FILE__ ) || isinteractive()
    test_reciprocity()
end