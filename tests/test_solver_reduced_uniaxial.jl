include("testconfig.jl")
using LinearAlgebra
using Statistics
using Test

@testset "Reduced Uniaxial Solver Tests" begin
    # Common parameters
    θ0 = 20.0  # incident angle in degrees
    λ = 632.8  # wavelength in nm
    
    # Surface parameters
    Nx = 2^12  # Number of discretization points
    Lx = 100 * λ  # Length of the surface
    
    # Material parameters
    ε_perp = -10.0 + 0.0im  # perpendicular permittivity (negative for high reflectivity)
    ε_para = -10.0 + 0.0im  # parallel permittivity
    μ_perp = 1.0 + 0.0im    # perpendicular permeability
    μ_para = 1.0 + 0.0im    # parallel permeability
    
    # Media
    above = Vacuum()

    @testset "Flat Surface - Uniaxial Medium" begin
        # Create flat surface and uniaxial medium
        surface = FlatSurface()
        below = Uniaxial(ε_perp, ε_para, μ_perp, μ_para)
        
        # Create parameters and solver data
        params = Parameters(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=[θ0])
        
        # Create solver data
        data = SolverData(params, 1, :reduced)
        
        # Preallocate and precompute
        precomputed = Precomputed(data)
        precompute!(precomputed, data)
        
        # Run the solver for a single realization
        alloc = Preallocated(data)
        generate_surface!(alloc, params)
        solve_single!(alloc, precomputed, data)
        
        # Record results for observation
        observe!(data, alloc, 1)
        
        # For a flat surface with negative permittivity, we should have high reflectivity
        k_idx = 1  # Only one incident angle
        
        # Find the index closest to the specular angle
        sin_theta = sind(θ0)
        _, q_idx = findmin(abs.(params.qs .- sin_theta))
        
        # Verify reflection coefficients are high for this reflective material
        r_p = abs(alloc.PNpk[q_idx, k_idx])
        r_s = abs(alloc.SNpk[q_idx, k_idx])
        
        @test r_p > 0.9  # Expect high reflection for reflective material
        @test r_s > 0.9  # Expect high reflection for reflective material
        
        # Check energy conservation
        P_energy, S_energy = energy_conservation(data)
        @test isapprox(P_energy[1], 1.0, rtol=1e-2)  # Energy should be conserved
        @test isapprox(S_energy[1], 1.0, rtol=1e-2)  # Energy should be conserved
        
        # For normal incidence, the p and s polarizations should have the same reflection coefficient
        if isapprox(θ0, 0.0, atol=1e-6)
            k_idx = 1  # Only one incident angle
            _, q_idx = findmin(abs.(params.qs .- 0.0))
            
            r_p = alloc.PNpk[q_idx, k_idx]
            r_s = alloc.SNpk[q_idx, k_idx]
            
            @test isapprox(abs(r_p), abs(r_s), rtol=1e-3)
        end
    end
    
    # Test random Gaussian surface with uniaxial medium
    @testset "Gaussian Surface - Uniaxial Medium" begin
        # Create a Gaussian surface with small roughness
        σ = 10.0e-9  # RMS height (smaller for better numerical stability)
        ℓ = 100.0e-9  # Correlation length
        surface = GaussianSurface(σ, ℓ)
        
        # Use the same uniaxial material
        below = Uniaxial(ε_perp, ε_para, μ_perp, μ_para)
        
        # More incident angles for a more comprehensive test
        θs = [10.0, 20.0, 30.0, 40.0, 50.0]
        
        params = Parameters(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)
        
        # Create solver data
        iters = 5  # Small number for test purposes
        data = SolverData(params, iters, :reduced)
        
        # Solve the ensemble problem
        solve_ensemble!(data)
        
        # Calculate MDRC
        mdrc = calc_mdrc(data)
        
        # Check that all reflectances are non-negative
        @test all(mdrc.coh_p .>= 0.0)
        @test all(mdrc.coh_s .>= 0.0)
        @test all(mdrc.inc_p .>= 0.0)
        @test all(mdrc.inc_s .>= 0.0)
    end
    
    # Test for reflection with varying anisotropy ratios
    @testset "Reflection - Varying Anisotropy" begin
        σ = 10.0e-9
        ℓ = 100.0e-9
        surface = GaussianSurface(σ, ℓ)
        
        # Test different anisotropy ratios using reflective materials
        anisotropy_ratios = [1.0, 1.5, 2.0, 3.0]
        
        for ratio in anisotropy_ratios
            # Create uniaxial medium with specified anisotropy ratio
            # Using the base ε_perp but varying the ratio for ε_para
            anisotropic_para = ratio * ε_perp
            
            below = Uniaxial(ε_perp, anisotropic_para, μ_perp, μ_para)
            
            params = Parameters(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=[θ0])
            
            # Create solver data
            data = SolverData(params, 1, :reduced)
            
            # Preallocate and precompute
            precomputed = Precomputed(data)
            precompute!(precomputed, data)
            
            # Run the solver for a single realization
            alloc = Preallocated(data)
            generate_surface!(alloc, params)
            solve_single!(alloc, precomputed, data)
            
            # Record results for observation
            observe!(data, alloc, 1)
            
            # Check energy conservation
            P_energy, S_energy = energy_conservation(data)
            @test isapprox(P_energy[1], 1.0, rtol=1e-2)
            @test isapprox(S_energy[1], 1.0, rtol=1e-2)
            @info "Anisotropy: $(ratio), P_energy: $(P_energy[1]), S_energy: $(S_energy[1])"
        end
    end
    
    # Test that μ_perp = μ_para assertion is enforced in the reduced solver
    @testset "Assertion for μ_perp = μ_para" begin
        # Create a Gaussian surface
        σ = 10.0e-9
        ℓ = 100.0e-9
        surface = GaussianSurface(σ, ℓ)
        
        # Create uniaxial medium with unequal permeability
        different_mu = 1.2 + 0.0im      # Different mu parallel
        
        below = Uniaxial(ε_perp, ε_para, μ_perp, different_mu)
        
        # Create parameters with solver_method = :reduced
        # The assertion should fail during validation
        params = Parameters(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=[θ0])
        
        # When creating Parameters with solver_method=:reduced and Uniaxial with μ_perp ≠ μ_para,
        # the assertion should fail during precomputation
        data = SolverData(params, 1, :reduced)
        precomputed = Precomputed(data)
        @test_throws AssertionError precompute!(precomputed, data)
    end
end
