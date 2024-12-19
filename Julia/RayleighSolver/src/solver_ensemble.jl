function solve_MDRC!(data::SolverData{_P}) where {_P}
    # Solve the Rayleigh problem for a given Parameters struct
    # sp is the preallocated surface struct, containing the preallocated output, this is mutated in place
    # params is the Parameters struct, containing the constant parameters for the calculation
    # surf_generator! is a function that generates the surface sp.ys for each iteration
    # N_ens is the number of ensemble averages to perform
    # returns the coherent and incoherent MDRC

    @info "precompute!:"
    @time precompute!(data.precomputed, data.params)
    validate(data.precomputed)

    Nthreads = Threads.nthreads()
    LinearAlgebra.BLAS.set_num_threads(1)
    Niters = data.iters

    enable_threads = true

    @time begin

    if Niters > Nthreads && enable_threads
        @info "Mainloop: Running $Niters iters on $Nthreads threads..."
        
        allocs = [Preallocated(data.params) for _ in 1:Nthreads]
        lk = ReentrantLock()

        Threads.@threads :static for n in 1:Niters
            tidx = Threads.threadid()
            _alloc = allocs[tidx]

            generate_surface!(_alloc, data.params)
            solve_single!(_alloc, data)
            @lock lk begin
                observe!(data.P_res, _alloc.PNpk, n)
                observe!(data.S_res, _alloc.SNpk, n)
            end # lock
        end
    else
        @info "Mainloop: Running $Niters iters sequentially."
        alloc = Preallocated(data.params)
        for n in 1:Niters
            generate_surface!(alloc, data.params)
            solve_single!(alloc, data)
            observe!(alloc, n)
        end
        combine!(data.P_res, alloc.P_res, 1, 1)
        combine!(data.S_res, alloc.S_res, 1, 1)
    end

    end # @time


    return nothing
end

function get_coh_inc(R, R2, qs::Vector{Float64}, ks::Vector{Float64})
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

struct DataMDRC
    coh_p::Matrix{Float64}
    inc_p::Matrix{Float64}
    coh_s::Matrix{Float64}
    inc_s::Matrix{Float64}
    θis::Vector{Float64}
    θss::Vector{Float64}
end

"Returns coherent and incoherent MDRC calculated from ⟨R⟩ and ⟨R²⟩"
function calc_mdrc(data::SolverData)
    params = data.params
    qs = params.qs
    ks = params.ks

    mask = qs .>= -1.0 .&& qs .<= 1.0
    θss = asind.(qs[mask])
    θis = data.params.θs

    @assert all(ks .≈ sind.(θis))
    Rp = data.P_res.R
    R2p = data.P_res.R²

    Rs = data.S_res.R
    R2s = data.S_res.R²

    coh_p, inc_p = get_coh_inc(Rp[mask, :], R2p[mask, :], qs[mask], ks)
    coh_s, inc_s = get_coh_inc(Rs[mask, :], R2s[mask, :], qs[mask], ks)

    return DataMDRC(coh_p, inc_p, coh_s, inc_s, θis, θss)
end