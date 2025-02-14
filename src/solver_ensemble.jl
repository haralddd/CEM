function assert_branch_cuts(qs, material)
    for q in qs
        @assert imag(alpha(q, material)) >= 0.0
        @assert real(alpha(q, material)) >= 0.0
    end
end



function solve_MDRC!(data::SolverData{_P}) where {_P}
    # Solve the Rayleigh problem for a given Parameters struct
    # sp is the preallocated surface struct, containing the preallocated output, this is mutated in place
    # params is the Parameters struct, containing the constant parameters for the calculation
    # surf_generator! is a function that generates the surface sp.ys for each iteration
    # N_ens is the number of ensemble averages to perform
    # returns the coherent and incoherent MDRC

    @info "Running tests:"
    assert_branch_cuts(data.params.qs, data.params.above)
    assert_branch_cuts(data.params.qs, data.params.below)
    @info "Branch cuts OK"

    @info "Precomputing matrix elements:"
    @time precompute!(data.precomputed, data.params)
    validate(data.precomputed)
    do_debug = get(ENV, "JULIA_DEBUG", "") != ""
    show_iter = get(ENV, "JULIA_SHOWITERS", "") == "true"
    @info "debug: $do_debug"
    @info "show iters: $show_iter"

    Nthreads = Threads.nthreads()
    LinearAlgebra.BLAS.set_num_threads(1)
    Niters = data.iters

    enable_threads = Nthreads > 1
    iter = 0

    if Niters > Nthreads && enable_threads
        @info "Mainloop: Running $Niters iters on $Nthreads threads..."
        
        allocs = [Preallocated(data.params) for _ in 1:Nthreads]
        lk = ReentrantLock()

        @time Threads.@threads :static for _ in 1:Niters
            tidx = Threads.threadid()
            _alloc = allocs[tidx]

            generate_surface!(_alloc, data.params)
            solve_single!(_alloc, data)
            @lock lk begin
                iter += 1
                observe!(data.P_res, _alloc.PNpk, iter)
                observe!(data.S_res, _alloc.SNpk, iter)
                if show_iter
                    @info "Iter: $iter / $Niters"
                end
                if do_debug
                    P_energy, S_energy = energy_conservation(abs2(_alloc.PNpk), abs2(_alloc.SNpk), data.params)
                    @debug "P_energy: $P_energy"
                    @debug "S_energy: $S_energy"
                end
            end # lock
        end
    else
        @info "Mainloop: Running $Niters iters sequentially."
        alloc = Preallocated(data.params)
        @time for n in 1:Niters
            generate_surface!(alloc, data.params)
            solve_single!(alloc, data)
            observe!(data.P_res, alloc.PNpk, n)
            observe!(data.S_res, alloc.SNpk, n)
            iter += 1
            if show_iter
                @info "Iter: $iter / $Niters"
            end
            if do_debug
                P_energy, S_energy = energy_conservation(abs2(alloc.PNpk), abs2(alloc.SNpk), data.params)
                @debug "P_energy: $P_energy"
                @debug "S_energy: $S_energy"
            end
        end
    end


    return nothing
end

function get_mdrc_coh_inc(R, R2, qs::Vector{Float64}, ks::Vector{Float64}, params::Parameters)
    Lx = params.Lx
    λ = params.lambda

    coh = Matrix{Float64}(undef, (length(qs), length(ks)))
    inc = Matrix{Float64}(undef, (length(qs), length(ks)))
    for (j, k) in enumerate(ks)
        for (i, q) in enumerate(qs)
            C = abs2(alpha0(q)) / abs(alpha0(k))
            coh[i, j] = C * abs2(R[i, j])
            inc[i, j] = C * R2[i, j] - coh[i, j]
        end
    end
    return coh, inc
end

function get_mdtc_coh_inc(T, T2, qs::Vector{Float64}, ks::Vector{Float64}, params::Parameters, ν::Symbol = :p)
    Lx = params.Lx
    λ = params.lambda

    coh = Matrix{Float64}(undef, (length(qs), length(ks)))
    inc = Matrix{Float64}(undef, (length(qs), length(ks)))
    κpa = ν == :p ? params.below.eps_para : params.below.mu_para
    for (j, k) in enumerate(ks)
        for (i, q) in enumerate(qs)
            C = abs(alpha(q, params.below)) * abs(alpha0(q)) / (κpa * abs(alpha0(k)))
            coh[i, j] = C * abs2(T[i, j])
            inc[i, j] = C * T2[i, j] - coh[i, j]
        end
    end
    return coh, inc
end

struct PlotData
    coh_p::Matrix{Float64}
    inc_p::Matrix{Float64}
    coh_s::Matrix{Float64}
    inc_s::Matrix{Float64}
    θs::Vector{Float64}
    θ0s::Vector{Union{Float64,Tuple{Float64,Float64}}}
end

"Returns coherent and incoherent MDRC calculated from ⟨R⟩ and ⟨R²⟩"
function calc_mdrc(data::SolverData)
    params = data.params
    qs = params.qs
    ks = params.ks

    mask = qs .>= -1.0 .&& qs .<= 1.0
    θss = asind.(qs[mask])
    θ0s = data.params.θs

    @assert all(ks .≈ sind.(θis))
    Rp = data.P_res.R
    R2p = data.P_res.R²

    Rs = data.S_res.R
    R2s = data.S_res.R²
    Lx = data.params.Lx

    coh_p, inc_p = get_mdrc_coh_inc(Rp[mask, :], R2p[mask, :], qs[mask], ks, Lx)
    coh_s, inc_s = get_mdrc_coh_inc(Rs[mask, :], R2s[mask, :], qs[mask], ks, Lx)

    return PlotData(coh_p, inc_p, coh_s, inc_s, θss, θ0s)
end

"Returns coherent and incoherent MDRC calculated from ⟨R⟩ and ⟨R²⟩"
function calc_mdrc(data::SolverData{Parameters{_S,Vacuum,Uniaxial}}) where _S
    params = data.params
    qs = params.qs
    ks = params.ks

    mask = qs .>= -1.0 .&& qs .<= 1.0
    θss = asind.(qs[mask])
    θ0s = data.params.θs

    @assert all(ks .≈ sind.(θ0s))
    Rp = get_R(data.P_res)
    R2p = get_R²(data.P_res)

    Rs = get_R(data.S_res)
    R2s = get_R²(data.S_res)

    coh_p, inc_p = get_mdrc_coh_inc(Rp[mask, :], R2p[mask, :], qs[mask], ks, params)
    coh_s, inc_s = get_mdrc_coh_inc(Rs[mask, :], R2s[mask, :], qs[mask], ks, params)

    dq = params.dq
    @info "∑MDRC_s = $((sum(coh_s) + sum(inc_s))*dq)"
    @info "∑MDRC_p = $((sum(coh_p) + sum(inc_p))*dq)"

    return PlotData(coh_p, inc_p, coh_s, inc_s, θss, θ0s)
end



function θto(uni::Uniaxial, θ0::Float64)
    μpa = uni.mu_para
    εpa = uni.eps_para
    no = √(μpa*εpa)

    return asind(real(no)*sind(θ0))

end

function θte(uni::Uniaxial, θ0::Float64)
    εpe = uni.eps_perp
    εpa = uni.eps_para
    μpe = uni.mu_perp
    μpa = uni.mu_para

    ne2 = real(μpe*εpe)
    no2 = real(μpa*εpa)
    n(θ) = 1/√((cosd(θ)^2)/no2 + (sind(θ)^2)/ne2)
    f(θ) = n(θ)*sind(θ) - sind(θ0)

    return find_zero(f, θ0)
end

"Returns coherent and incoherent MDTC calculated from ⟨T⟩ and ⟨T²⟩"
function calc_mdtc(data::SolverData{Parameters{_S,Vacuum,Uniaxial}}) where _S
    params = data.params
    qs = params.qs
    ks = params.ks
    below = data.params.below

    mask = qs .>= -1.0 .&& qs .<= 1.0
    θts = asind.(qs[mask])
    θ0s = [(θto(below, θ0), θte(below, θ0)) for θ0 in params.θs]

    Tp = get_T(data.P_res)
    T2p = get_T²(data.P_res)

    Ts = get_T(data.S_res)
    T2s = get_T²(data.S_res)


    coh_p, inc_p = get_mdtc_coh_inc(Tp[mask, :], T2p[mask, :], qs[mask], ks, params, :p)
    coh_s, inc_s = get_mdtc_coh_inc(Ts[mask, :], T2s[mask, :], qs[mask], ks, params, :s)

    dq = params.dq
    @info "∑MDTC_s = $((sum(coh_s) + sum(inc_s))*dq)"
    @info "∑MDTC_p = $((sum(coh_p) + sum(inc_p))*dq)"

    return PlotData(coh_p, inc_p, coh_s, inc_s, θts, θ0s)
end