include("testconfig.jl")

using BenchmarkTools


function profile_gaussian_surfacegen()
    spa, sp = default_params_for_surface_testing(GaussianSurface(30.0e-9, 100.0e-9))
    function loop_de_loop(sp, spa)
        for _ in 1:1000
            generate_surface!(sp, spa)
        end
    end
    @btime generate_surface!($sp, $spa)
end

function profile_rectangular_surfacegen()
    spa, sp = default_params_for_surface_testing(RectangularSurface(30.0e-9, 0.82, 1.97))
    function loop_de_loop(sp, spa)
        for _ in 1:1000
            generate_surface!(sp, spa)
        end
    end
    @btime generate_surface!($sp, $spa)
end

function profile_isotropic_solver()

    @info "Simple glass isotropic:"
    data = config_glass_isotropic()
    @info "precompute"
    pc_stats = @benchmark precompute!($data)
    @info "Surface generation"
    surf_stats = @benchmark generate_surface!($data.sp, $data.spa)
    @info "Solve single"
    solve_single_stats = @benchmark solve_single!($data)
    @info "Observation"
    obs_stats = @benchmark observe!($data.out_p, $data.sp.p_data.Npk, 10)


    @info "Precomputation"
    display(pc_stats)
    @info "Surface generation"
    display(surf_stats)
    @info "Solve single"
    display(solve_single_stats)
    @info "Observation"
    display(obs_stats)
    
    # @code_warntype precompute!(data)
    # @code_warntype generate_surface!(data.sp, data.spa)
    # @code_warntype solve_single!(data)
    # @code_warntype observe!(data.out_p, data.sp.p_data.Npk, 11)
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
