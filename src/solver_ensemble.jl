"""
    assert_branch_cuts(qs, material)

Validates that the branch cuts are correctly defined for all wavenumbers in `qs`.
Ensures that the imaginary part of alpha is non-negative and the real part is non-negative.

# Arguments:
- `qs`: Vector of wavenumbers to check
- `material`: Material object containing the properties to calculate alpha
"""
function assert_branch_cuts(qs, material)
    for q in qs
        @assert imag(alpha(q, material)) >= 0.0
        @assert real(alpha(q, material)) >= 0.0
    end
end

function solve_single!(alloc::Preallocated, pre::Precomputed, data::SolverData)
    if data.solver_type == :full
        return solve_single_full!(alloc, pre, data)
    elseif data.solver_type == :reduced || data.solver_type == :combined
        return solve_single_reduced!(alloc, pre, data)
    end
end

"""
    solve_ensemble!(data::SolverData{_P}) where {_P}

Solves the Rayleigh problem for multiple surface realizations and computes ensemble averages.
Handles both multi-threaded and single-threaded computation depending on the system configuration.

# Arguments:
- `data`: [`SolverData`](@ref) - Contains the parameters, preallocated steps, and output structures
"""
function solve_ensemble!(data::SolverData{_P}) where {_P}
    @info "Precomputing matrix elements:"
    precomputed = Precomputed(data.params)
    @time precompute!(precomputed, data.params)
    validate(precomputed)
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
        
        allocs = [Preallocated(data) for _ in 1:Nthreads]
        lk = ReentrantLock()

        @time Threads.@threads :static for _ in (show_iter ? ProgressBar(1:Niters) : 1:Niters)
            tidx = Threads.threadid()
            _alloc = allocs[tidx]

            generate_surface!(_alloc, data.params)
            solve_single!(_alloc, precomputed, data)
            @lock lk begin
                iter += 1
                observe!(data, _alloc, iter)
            end # lock
        end
    else
        @info "Mainloop: Running $Niters iters sequentially."
        alloc = Preallocated(data)
        @time for n in (show_iter ? ProgressBar(1:Niters) : 1:Niters)
            generate_surface!(alloc, data.params)
            solve_single!(alloc, precomputed, data)
            observe!(data, alloc, n)
            iter += 1
        end
    end


    return nothing
end