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

tio2 = Uniaxial(6.84+0.01im, 3.62+0.01im, 1.0+0.0im, 1.0+0.0im)
function profile_ensemble_reduced()
    @info "Default uniaxial reduced ensemble solve"
    paramconf = ParametersConfig(
        Nx=4096,
        Lx=200,
        below=tio2,
    )
    data = SolverData(paramconf, 1000, :reduced)
    @btime solve_ensemble!($data)
end

function profile_ensemble_full()
    @info "Default uniaxial full ensemble solve"
    paramconf = ParametersConfig(
        Nx=4096,
        Lx=200,
        below=tio2,
    )
    data = SolverData(paramconf, 1000, :full)
    @btime solve_ensemble!($data)
end
    
