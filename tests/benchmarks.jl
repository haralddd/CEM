include("testconfig.jl")

using BenchmarkTools
# using ProfileView

function profile_gaussian_surfacegen()
    @info "Gaussian surface generation"
    data = default_gaussian_config()
    pre = Preallocated(data.params)
    bench = @benchmark generate_surface!($pre, $data.params)
    display(bench)
end

function profile_rectangular_surfacegen()
    @info "Rectangular surface generation"
    data = default_rectangular_config()
    pre = Preallocated(data.params)
    
    bench = @benchmark generate_surface!($pre, $data.params)
    display(bench)
end

function profile_isotropic_precompute()
    @info "Simple glass isotropic precompute"
    data = config_glass_isotropic()

    bench = @benchmark precompute!($data)
    display(bench)
end

function profile_isotropic_single_solve()
    @info "Simple glass isotropic single solve"
    data = config_glass_isotropic()
    precompute!(data)
    generate_surface!(data.sp, data.params)

    bench = @benchmark solve_single!($data)
    display(bench)
end

function profile_isotropic_observe()
    @info "Simple glass isotropic observe"
    data = config_glass_isotropic()
    precompute!(data)
    generate_surface!(data.sp, data.params)
    solve_single!(data)
    bench = @benchmark observe!($data, 10)
    display(bench)
end

function singlerun_isotropic_single_solve()
    data = config_glass_isotropic()
    precompute!(data)
    generate_surface!(data.sp, data.params)
    _, stats... = @timed solve_single!(data)
    @info "Single Solve Total:"
    display(stats.time)
end