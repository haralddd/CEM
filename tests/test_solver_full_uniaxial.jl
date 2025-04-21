include("testconfig.jl")
using CairoMakie
using LinearAlgebra
using Statistics
using Test

@testset "Full Uniaxial Solver Tests" begin
    # Common parameters
    θ0 = 20.0  # incident angle in degrees
    λ = 632.8  # wavelength in nm
    
    # Surface parameters
    Nx = 2^10  # Number of discretization points
    Lx = 50 * λ  # Length of the surface
    
    # Media
    above = Vacuum()

    @testset "Flat Surface - Isotropic Uniaxial Medium" begin
        # Create uniaxial medium with isotropic properties (εₚₑᵣₚ = εₚₐᵣₐ, μₚₑᵣₚ = μₚₐᵣₐ)
        ε_perp = 2.0 + 0.0im
        ε_para = 2.0 + 0.0im
        μ_perp = 1.0 + 0.0im
        μ_para = 1.0 + 0.0im
        
        # Create flat surface and uniaxial medium
        surface = FlatSurface()
        below = Uniaxial(ε_perp, ε_para, μ_perp, μ_para)
        
        # Create parameters and solver data - using full solver
        params = Parameters(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=[θ0])
        data = SolverData(params, 1, :full)
        
        # Preallocate and precompute
        precomputed = Precomputed(params)
        precompute_full!(precomputed, params)
        validate(precomputed)
        
        # Run the solver for a single realization
        alloc = Preallocated(data)
        generate_surface!(alloc, params)
        solve_single!(alloc, precomputed, data)
        observe!(data, alloc, 1)
        
        # For a flat surface with isotropic properties, we should have predictable reflection/transmission
        # Find the specular reflection/transmission index
        k_idx = 1  # Only one incident angle
        
        # Find the index closest to the specular angle
        sin_theta = sind(θ0)
        _, q_idx = findmin(abs.(params.qs .- sin_theta))
        
        # For an isotropic medium, we can verify against Fresnel's equations
        # Calculate expected reflection/transmission coefficients
        n1 = 1.0  # Vacuum
        n2 = sqrt(ε_perp * μ_perp)  # Isotropic uniaxial medium
        θi = deg2rad(θ0)
        
        # Calculate refraction angle using Snell's law
        θt = asin(sin(θi) / real(n2))
        
        # Calculate Fresnel coefficients for p-polarization (TM)
        r_p_expected = (n2 * cos(θi) - n1 * cos(θt)) / (n2 * cos(θi) + n1 * cos(θt))
        t_p_expected = (2 * n1 * cos(θi)) / (n2 * cos(θi) + n1 * cos(θt))
        
        # Calculate Fresnel coefficients for s-polarization (TE)
        r_s_expected = (n1 * cos(θi) - n2 * cos(θt)) / (n1 * cos(θi) + n2 * cos(θt))
        t_s_expected = (2 * n1 * cos(θi)) / (n1 * cos(θi) + n2 * cos(θt))
        
        # Extract actual coefficients from solver results
        # For full solver, the first half of the alloc.PNpk is reflection, second half is transmission
        half_size = size(alloc.PNpk, 1) ÷ 2
        
        # Reflection coefficients
        r_p_actual = alloc.PNpk[q_idx, k_idx]
        r_s_actual = alloc.SNpk[q_idx, k_idx]
        
        # Transmission coefficients
        t_p_actual = alloc.PNpk[half_size + q_idx, k_idx]
        t_s_actual = alloc.SNpk[half_size + q_idx, k_idx]
        
        # Test that calculated values match expected values
        @test isapprox(abs(r_p_actual), abs(r_p_expected), rtol=1e-2)
        @test isapprox(abs(r_s_actual), abs(r_s_expected), rtol=1e-2)
        @test isapprox(abs(t_p_actual), abs(t_p_expected), rtol=1e-2)
        @test isapprox(abs(t_s_actual), abs(t_s_expected), rtol=1e-2)
        
        # Check energy conservation
        P_energy, S_energy = energy_conservation(data)
        @test isapprox(P_energy[1], 1.0, rtol=1e-2)
        @test isapprox(S_energy[1], 1.0, rtol=1e-2)
        @info "P_energy: $(P_energy[1]), S_energy: $(S_energy[1])"
    end
    
    @testset "Flat Surface - Anisotropic Uniaxial Medium" begin
        # Create uniaxial medium with anisotropic permittivity (εₚₑᵣₚ ≠ εₚₐᵣₐ)
        ε_perp = 2.0 + 0.0im
        ε_para = 3.0 + 0.0im  # Different from ε_perp
        μ_perp = 1.0 + 0.0im
        μ_para = 1.0 + 0.0im
        
        # Create flat surface and uniaxial medium
        surface = FlatSurface()
        below = Uniaxial(ε_perp, ε_para, μ_perp, μ_para)
        
        # Create parameters and solver data
        params = Parameters(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=[θ0])
        data = SolverData(params, 1, :full)
        
        # Preallocate and precompute
        precomputed = Precomputed(params)
        precompute_full!(precomputed, params)
        validate(precomputed)
        
        # Run the solver for a single realization
        alloc = Preallocated(data)
        generate_surface!(alloc, params)
        solve_single!(alloc, precomputed, data)
        observe!(data, alloc, 1)
        
        # Check energy conservation for anisotropic medium
        P_energy, S_energy = energy_conservation(data)
        @test isapprox(P_energy[1], 1.0, rtol=1e-2)
        @test isapprox(S_energy[1], 1.0, rtol=1e-2)
        @info "P_energy: $(P_energy[1]), S_energy: $(S_energy[1])"
    end
    
    @testset "Flat Surface - Anisotropic Permeability" begin
        # Create uniaxial medium with anisotropic permeability (μₚₑᵣₚ ≠ μₚₐᵣₐ)
        # This is valid for the full solver but not for the reduced solver
        ε_perp = 2.0 + 0.0im
        ε_para = 2.0 + 0.0im
        μ_perp = 1.0 + 0.0im
        μ_para = 1.2 + 0.0im  # Different from μ_perp
        
        # Create flat surface and uniaxial medium
        surface = FlatSurface()
        below = Uniaxial(ε_perp, ε_para, μ_perp, μ_para)
        
        # Create parameters and solver data
        params = Parameters(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=[θ0])
        data = SolverData(params, 1, :full)
        
        # Preallocate and precompute
        precomputed = Precomputed(params)
        precompute_full!(precomputed, params)
        validate(precomputed)
        
        # Run the solver for a single realization
        alloc = Preallocated(data)
        generate_surface!(alloc, params)
        solve_single!(alloc, precomputed, data)
        observe!(data, alloc, 1)
        
        # Check energy conservation with anisotropic permeability
        P_energy, S_energy = energy_conservation(data)
        @test isapprox(P_energy[1], 1.0, rtol=1e-2)
        @test isapprox(S_energy[1], 1.0, rtol=1e-2)
        @info "P_energy: $(P_energy[1]), S_energy: $(S_energy[1])"
    end
    
    @testset "Gaussian Surface - Uniaxial Medium" begin
        # Create a Gaussian surface with small roughness
        σ = 0.05 * λ  # RMS height (smaller for better numerical stability)
        ℓ = 1.0 * λ  # Correlation length
        surface = GaussianSurface(σ, ℓ)
        
        # Create uniaxial medium with anisotropic properties
        ε_perp = 2.0 + 0.0im
        ε_para = 3.0 + 0.0im
        μ_perp = 1.0 + 0.0im
        μ_para = 1.2 + 0.0im
        below = Uniaxial(ε_perp, ε_para, μ_perp, μ_para)
        
        # More incident angles for a more comprehensive test
        θs = [10.0, 20.0, 30.0, 40.0, 50.0]
        params = Parameters(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)
        
        # Create solver data
        iters = 5  # Small number for test purposes
        data = SolverData(params, iters, :full)
        
        # Solve the ensemble problem
        solve_ensemble!(data)
        
        # Calculate MDRC and MDTC
        mdrc = calc_mdrc(data)
        mdtc = calc_mdtc(data)
        
        # Check that all reflectances and transmittances are non-negative
        @test all(mdrc.coh_p .>= 0.0)
        @test all(mdrc.coh_s .>= 0.0)
        @test all(mdrc.inc_p .>= 0.0)
        @test all(mdrc.inc_s .>= 0.0)
        
        @test all(mdtc.coh_p .>= 0.0)
        @test all(mdtc.coh_s .>= 0.0)
        @test all(mdtc.inc_p .>= 0.0)
        @test all(mdtc.inc_s .>= 0.0)
        
        # Check energy conservation
        P_energy, S_energy = energy_conservation(data)
        for i in eachindex(P_energy)
            @test isapprox(P_energy[i], 1.0, rtol=1e-2)
            @test isapprox(S_energy[i], 1.0, rtol=1e-2)
        end
        @info "P_energy: $(P_energy), S_energy: $(S_energy)"
    end
end
