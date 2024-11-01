include("testconfig.jl")

using ProfileView
using BenchmarkTools


function profile_gaussian_surfacegen()
    spa, sp = default_params_for_surface_testing(GaussianSurface(30.0e-9, 100.0e-9))
    function loop_de_loop(sp, spa)
        for _ in 1:1000
            generate_surface!(sp, spa)
        end
    end
    @btime generate_surface!($sp, $spa)
    ProfileView.@profview loop_de_loop(sp, spa)
end

function profile_rectangular_surfacegen()
    spa, sp = default_params_for_surface_testing(RectangularSurface(30.0e-9, 0.82, 1.97))
    function loop_de_loop(sp, spa)
        for _ in 1:1000
            generate_surface!(sp, spa)
        end
    end
    @btime generate_surface!($sp, $spa)
    ProfileView.@profview loop_de_loop(sp, spa)
end

function profile_isotropic_solver()

    @info "Simple glass isotropic:"
    spa = config_glass_isotropic()
    data = SolverData(spa, 100)
    _, (pc_stats...) = @timed precompute!(data)
    _, (surf_stats...) = @timed generate_surface!(data.sp, data.spa)
    _, (solve_single_stats...) = @timed solve_single!(data)
    _, (obs_stats...) = @timed observe!(data.out, data.sp.Npk, 10)

    @info "Precomputation: $pc_stats"
    @info "Surface generation: $surf_stats"
    @info "Single solve: $solve_single_stats"
    @info "Observation: $obs_stats"

    @info "Simple silver isotropic"
    spa = config_silver_isotropic()
    data = SolverData(spa, 100)
    _, (pc_stats...) = @timed precompute!(data)
    _, (surf_stats...) = @timed generate_surface!(data.sp, spa)
    _, (solve_single_stats...) = @timed solve_single!(data)
    _, (obs_stats...) = @timed observe!(data.out, data.sp.Npk, 10)

    @info "Precomputation: $pc_stats"
    @info "Surface generation: $surf_stats"
    @info "Single solve: $solve_single_stats"
    @info "Observation: $obs_stats"
end

profile_isotropic_solver()

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
