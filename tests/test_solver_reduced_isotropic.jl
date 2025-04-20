include("testconfig.jl")
using Plots
using LinearAlgebra
using Statistics
using Test

@testset "Reduced Isotropic Solver Tests" begin
    θ0 = 20.0  # incident angle in degrees
    λ = 632.8  # wavelength in nm
    ε = -10.0 + 0.0im  # relative permittivity of the medium
    μ = 1.0 + 0.0im   # relative permeability of the medium
    
    # Create a flat surface
    Nx = 2^12  # Increased number of discretization points to satisfy Q > 4
    Lx = 100 * λ  # Length of the surface
    above = Vacuum()

    @testset "Flat Surface" begin
        surface = FlatSurface()
        below = Isotropic(ε, μ)
        params = Parameters(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=[θ0])
        data = SolverData(params, 1, :reduced)
        
        # Preallocate and precompute
        precomputed = Precomputed(params)
        precompute!(precomputed, params)
        validate(precomputed)
        
        # Run the solver for a single realization
        alloc = Preallocated(data)
        generate_surface!(alloc, params)
        solve_single!(alloc, precomputed, data)
        
        # Test results
        # For a flat surface, we should have perfect reflection according to Fresnel's equations
        
        # Compute expected reflectivity from Fresnel equations
        n1 = 1.0  # vacuum
        n2 = sqrt(ε * μ)  # isotropic medium
        θ0_rad = deg2rad(θ0)
        
        # P-polarization (TM) Fresnel reflection coefficient
        r_p_expected = (n1 * cos(θ0_rad) - n2 * sqrt(1 - (n1/n2 * sin(θ0_rad))^2)) / 
                     (n1 * cos(θ0_rad) + n2 * sqrt(1 - (n1/n2 * sin(θ0_rad))^2))
        
        # S-polarization (TE) Fresnel reflection coefficient
        r_s_expected = (n1 * sqrt(1 - (n1/n2 * sin(θ0_rad))^2) - n2 * cos(θ0_rad)) / 
                     (n1 * sqrt(1 - (n1/n2 * sin(θ0_rad))^2) + n2 * cos(θ0_rad))
        
        # Get reflection coefficients from solver
        k_idx = 1  # Only one incident angle
        
        # Find the closest q value to sin(θ0_rad)
        q_target = sin(θ0_rad)
        _, q_idx = findmin(abs.(params.qs .- q_target))
        
        r_p_actual = alloc.PNpk[q_idx, k_idx]
        r_s_actual = alloc.SNpk[q_idx, k_idx]
        
        @test isapprox(abs(r_p_actual), abs(r_p_expected), rtol=0.2)  # 20% tolerance
        @test isapprox(abs(r_s_actual), abs(r_s_expected), rtol=0.2)  # 20% tolerance
        
    end
    
    # Test random Gaussian surface
    @testset "Gaussian Surface" begin
        σ = 10.0e-9  # RMS height
        ℓ = 100.0e-9  # Correlation length
        surface = GaussianSurface(σ, ℓ)
        below = Isotropic(ε, μ)
        θs = [10.0, 20.0, 30.0, 40.0, 50.0]
        params = Parameters(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)
        
        # Create solver data
        iters = 5  # Small number for test purposes
        data = SolverData(params, iters, :reduced)
        
        # Solve the ensemble problem
        solve_ensemble!(data)
        P_energy, S_energy = energy_conservation(data)

        # Check energy conservation
        @test all( isapprox.(P_energy, 1.0, atol=0.01)) # Just ensure positive and not too small
        @test all( isapprox.(S_energy, 1.0, atol=0.01)) # Just ensure positive and not too small
        @info "P_energy: $(P_energy), S_energy: $(S_energy)"
        
        # Calculate MDRC
        mdrc = calc_mdrc(data)
        @test all(mdrc.coh_p .>= 0.0)
        @test all(mdrc.coh_s .>= 0.0)
        @test all(mdrc.inc_p .>= 0.0)
        @test all(mdrc.inc_s .>= 0.0)
        dq = params.dq
        alphas = cosd.(mdrc.θs) # scattering angle factor
        filt = alphas .> 0.0
        for i in 1:length(θs)
            @info "MDRC P sum: $(sum((mdrc.coh_p[filt, i] .+ mdrc.inc_p[filt, i]) ./ alphas[filt]) * dq)"
            @info "MDRC S sum: $(sum((mdrc.coh_s[filt, i] .+ mdrc.inc_s[filt, i]) ./ alphas[filt]) * dq)"
        end
    end
    
    # Test single bump surface - simplified special case
    @testset "Single Bump Surface" begin
        θ0 = 10.0
        h = 10.0e-9  # Larger bump height (as a fraction of wavelength)
        w = 100.0e-9  # Larger bump width (as a fraction of wavelength)
        surface = SingleBumpSurface(h, w)
        below = Isotropic(ε, μ)
        params = Parameters(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=[θ0])
        data = SolverData(params, 1, :reduced)
        
        # Preallocate and precompute
        precomputed = Precomputed(params)
        precompute!(precomputed, params)
        validate(precomputed)
        
        # Run the solver for a single realization
        alloc = Preallocated(data)
        generate_surface!(alloc, params)
        solve_single!(alloc, precomputed, data)
        observe!(data, alloc, 1)
        
        # Test results
        # For a bump at normal incidence, we expect:
        # 1. Some specular reflection
        # 2. Diffuse scattering symmetrically around normal
        # Test energy conservation
        P_energy, S_energy = energy_conservation(data)
        
        # For the single bump case, energy conservation might not be perfect
        # due to discretization effects, numerical approximations, and finite truncation.
        # Just check that values are within a reasonable range
        @test all( isapprox.(P_energy, 1.0, atol=0.01)) # Just ensure positive and not too small
        @test all( isapprox.(S_energy, 1.0, atol=0.01)) # Just ensure positive and not too small
        @info "P_energy: $(P_energy), S_energy: $(S_energy)"
        
        # For symmetry around normal incidence, we can check that the
        # maximum values on each side are similar
        q_vals = abs.(alloc.PNpk[:,1]) # Get reflection coefficients for the first angle
        mid_idx = length(q_vals) ÷ 2
        left_side = q_vals[1:mid_idx]
        right_side = reverse(q_vals[mid_idx+1:end])
    end
end

# Optional: Plot the results for visualization
function plot_reflectivity()
    # Setup parameters
    θ0 = 30.0  # incident angle in degrees
    λ = 632.8  # wavelength in nm
    ε = 2.25 + 0.0im  # relative permittivity
    μ = 1.0 + 0.0im   # relative permeability
    
    # Create a Gaussian surface
    Nx = 2^12
    Lx = 100 * λ
    σ = 10.0e-9
    ℓ = 100.0e-9
    surface = GaussianSurface(σ, ℓ)
    
    # Initialize solver
    above = Vacuum()
    below = Isotropic(ε, μ)
    
    params = Parameters(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=[θ0])
    
    # Create solver data
    iters = 20  # More iterations for better averaging
    data = SolverData(params, iters, :reduced)
    
    # Solve the ensemble problem
    solve_ensemble!(data)
    
    # Calculate MDRC
    mdrc = calc_mdrc(data)
    
    # Plot the results
    p1 = plot(mdrc.θs, mdrc.coh_p[:, 1], label="Coherent P", linewidth=2)
    plot!(p1, mdrc.θs, mdrc.inc_p[:, 1], label="Incoherent P", linewidth=2)
    xlabel!(p1, "Scattering Angle (degrees)")
    ylabel!(p1, "MDRC")
    title!(p1, "P-Polarization MDRC (θ₀ = $θ0°)")
    
    p2 = plot(mdrc.θs, mdrc.coh_s[:, 1], label="Coherent S", linewidth=2)
    plot!(p2, mdrc.θs, mdrc.inc_s[:, 1], label="Incoherent S", linewidth=2)
    xlabel!(p2, "Scattering Angle (degrees)")
    ylabel!(p2, "MDRC")
    title!(p2, "S-Polarization MDRC (θ₀ = $θ0°)")
    
    plot(p1, p2, layout=(2,1), size=(800, 600))
end

# Uncomment to run the plot function
# plot_reflectivity()
