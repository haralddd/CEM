include("testconfig.jl")

using BenchmarkTools
# using ProfileView
function profile_uniaxial_reduced_single_solve()
    @info "Default uniaxial reduced single solve"
    data = config_default_uniaxial(solver_type=:reduced)
    pre = Precomputed(data)
    alloc = Preallocated(data)
    precompute!(pre, data)
    generate_surface!(alloc, data.params)

    @time solve_single!(alloc, pre, data)

    # bench = @benchmark solve_single_reduced!($alloc, $pre, $data)
    # display(bench)
end

function profile_uniaxial_hybrid_single_solve()
    @info "Default uniaxial hybrid single solve"
    data = config_default_uniaxial(solver_type=:hybrid)
    pre = Precomputed(data)
    alloc = Preallocated(data)
    precompute!(pre, data)
    generate_surface!(alloc, data.params)

    @time solve_single!(alloc, pre, data)

    # bench = @benchmark solve_single_hybrid!($alloc, $pre, $data)
    # display(bench)
end

function profile_uniaxial_full_single_solve()
    @info "Default uniaxial full single solve"
    data = config_default_uniaxial(solver_type=:full)
    pre = Precomputed(data)
    alloc = Preallocated(data)
    precompute!(pre, data)
    generate_surface!(alloc, data.params)

    @time solve_single!(alloc, pre, data)

    # bench = @benchmark solve_single!($data)
    # display(bench)
end