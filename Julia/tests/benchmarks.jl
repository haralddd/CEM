include("testconfig.jl")

using BenchmarkTools
using ProfileView

function profile_gaussian_surfacegen()
    @info "Gaussian surface generation"
    data = default_gaussian_config()
    bench = @benchmark generate_surface!($data.sp, $data.spa)
    display(bench)
end

function profile_rectangular_surfacegen()
    @info "Rectangular surface generation"
    data = default_rectangular_config()
    
    bench = @benchmark generate_surface!($data.sp, $data.spa)
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
    generate_surface!(data.sp, data.spa)

    bench = @benchmark solve_single!($data)
    display(bench)
end

function profile_crystal_solver()
    ε = 2.25 + 1e-4im
    lambda = 632.8e-9
    Q = 4
    Nq = 512 + 1
    ks = [sind(20.0)]
    L = 10.0e-6
    Ni = 3
    surf = GaussianSurface(30.0e-9, 100.0e-9)
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

    @benchmark SimPreCompute($rp_crystal)
end
