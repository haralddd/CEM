using RayleighSolver
using Test

@testset "Reduced Uniaxial Solver Tests" begin
    # Common parameters
    λ = 632.8  # wavelength in nm
    
    # Surface parameters
    Nx = 4096  # Number of discretization points
    Lx = 100  # Length of the surface
    
    # Material parameters
    ε_perp = -10.0 + 0.0im  # perpendicular permittivity (negative for high reflectivity)
    ε_para = -10.0 + 0.0im  # parallel permittivity
    μ_perp = 1.0 + 0.0im    # perpendicular permeability
    μ_para = 1.0 + 0.0im    # parallel permeability
    θs = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]
    
    # Media
    above = Vacuum()
    ensemble_iters = 1000

    @testset "Gaussian Surface - Isotropic Metal Energy" begin
        # Create Gaussian surface and uniaxial medium
        surface = GaussianSurface(10.0e-9, 100.0e-9)
        below = Uniaxial(ε_perp, ε_para, μ_perp, μ_para)
        
        # Create parameters and solver data
        params = ParametersConfig(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)
        
        # Create solver data
        data = SolverData(params, ensemble_iters, :reduced)
        solve_ensemble!(data)
        
        # Check energy conservation
        P_energy, S_energy = energy_conservation(data)
        @test isapprox(P_energy[1], 1.0, rtol=1e-2)  # Energy should be conserved
        @test isapprox(S_energy[1], 1.0, rtol=1e-2)  # Energy should be conserved
        @info "P_energy: $(P_energy[1]), S_energy: $(S_energy[1])"
    end
    
    # Test random Gaussian surface with uniaxial medium
    @testset "Gaussian Surface - Uniaxial Metal Energy 1" begin
        # Create a Gaussian surface with small roughness
        σ = 10.0e-9  # RMS height (smaller for better numerical stability)
        ℓ = 100.0e-9  # Correlation length
        surface = GaussianSurface(σ, ℓ)
        
        # Use the same uniaxial material
        below = Uniaxial(ε_perp, 0.5*ε_para, μ_perp, μ_para)
        
        params = ParametersConfig(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)
        
        # Create solver data
        data = SolverData(params, ensemble_iters, :reduced)
        
        # Solve the ensemble problem
        solve_ensemble!(data)
        
        # Calculate MDRC
        mdrc = calc_mdrc(data)
        
        # Check that all reflectances are non-negative
        @test all(mdrc.coh_p .>= 0.0)
        @test all(mdrc.coh_s .>= 0.0)
        @test all(mdrc.inc_p .>= 0.0)
        @test all(mdrc.inc_s .>= 0.0)

        P_energy, S_energy = energy_conservation(data)
        @test isapprox(P_energy[1], 1.0, rtol=1e-2)  # Energy should be conserved
        @test isapprox(S_energy[1], 1.0, rtol=1e-2)  # Energy should be conserved
        @info "P_energy: $(P_energy[1]), S_energy: $(S_energy[1])"
    end

        # Test random Gaussian surface with uniaxial medium
    @testset "Gaussian Surface - Uniaxial Metal Energy 1" begin
        # Create a Gaussian surface with small roughness
        σ = 10.0e-9  # RMS height (smaller for better numerical stability)
        ℓ = 100.0e-9  # Correlation length
        surface = GaussianSurface(σ, ℓ)
        
        # Use the same uniaxial material
        below = Uniaxial(ε_perp, 1.5*ε_para, μ_perp, μ_para)
        
        params = ParametersConfig(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)
        
        # Create solver data
        data = SolverData(params, ensemble_iters, :reduced)
        
        # Solve the ensemble problem
        solve_ensemble!(data)
        
        # Calculate MDRC
        mdrc = calc_mdrc(data)
        
        # Check that all reflectances are non-negative
        @test all(mdrc.coh_p .>= 0.0)
        @test all(mdrc.coh_s .>= 0.0)
        @test all(mdrc.inc_p .>= 0.0)
        @test all(mdrc.inc_s .>= 0.0)

        P_energy, S_energy = energy_conservation(data)
        @test isapprox(P_energy[1], 1.0, rtol=1e-2)  # Energy should be conserved
        @test isapprox(S_energy[1], 1.0, rtol=1e-2)  # Energy should be conserved
        @info "P_energy: $(P_energy[1]), S_energy: $(S_energy[1])"
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
        params = ParametersConfig(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)
        
        # When creating Parameters with solver_method=:reduced and Uniaxial with μ_perp ≠ μ_para,
        # the assertion should fail during precomputation
        data = SolverData(params, ensemble_iters, :reduced)
        @test_throws AssertionError solve_ensemble!(data)
    end

        @testset "Gaussian Surface - Bianisotropy HMM1" begin
        # Create a Gaussian surface with small roughness
        σ = 10.0e-9  # RMS height (smaller for better numerical stability)
        ℓ = 100.0e-9  # Correlation length
        surface = GaussianSurface(σ, ℓ)

        # Create uniaxial medium with anisotropic properties
        ε_perp = -7.0 + 0.0im
        ε_para = 3.0 + 0.0im
        μ_perp = 1.0 + 0.0im
        μ_para = 1.2 + 0.0im
        below = Uniaxial(ε_perp, ε_para, μ_perp, μ_para)

        # More incident angles for a more comprehensive test
        paramconf = ParametersConfig(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)

        # Create solver data
        data = SolverData(paramconf, ensemble_iters, :reduced)
        params = data.params

        # Solve the ensemble problem
        solve_ensemble!(data)

        # Calculate MDRC
        mdrc = calc_mdrc(data)

        # Check that all reflectances and transmittances are non-negative
        @test all(mdrc.coh_p .>= 0.0)
        @test all(mdrc.coh_s .>= 0.0)
        @test all(mdrc.inc_p .>= 0.0)
        @test all(mdrc.inc_s .>= 0.0)

        P_energy, S_energy = energy_conservation(data)
        for i in eachindex(P_energy)
            @test !isapprox(P_energy[i], 1.0, rtol=1e-2)  # Energy should not be conserved (transmission)
            @test !isapprox(S_energy[i], 1.0, rtol=1e-2)  # Energy should not be conserved (reflection)
        end
        @info "P_energy: $(P_energy), S_energy: $(S_energy)"
    end

    @testset "Gaussian Surface - Bianisotropy HMM2" begin
        # Create a Gaussian surface with small roughness
        σ = 10.0e-9  # RMS height (smaller for better numerical stability)
        ℓ = 100.0e-9  # Correlation length
        surface = GaussianSurface(σ, ℓ)

        # Create uniaxial medium with anisotropic properties
        ε_perp = 3.0 + 0.0im
        ε_para = -7.0 + 0.0im
        μ_perp = 1.0 + 0.0im
        μ_para = 1.2 + 0.0im
        below = Uniaxial(ε_perp, ε_para, μ_perp, μ_para)

        # More incident angles for a more comprehensive test
        paramconf = ParametersConfig(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)

        # Create solver data
        data = SolverData(paramconf, ensemble_iters, :reduced)
        params = data.params

        # Solve the ensemble problem
        solve_ensemble!(data)

        # Calculate MDRC
        mdrc = calc_mdrc(data)

        # Check that all reflectances and transmittances are non-negative
        @test all(mdrc.coh_p .>= 0.0)
        @test all(mdrc.coh_s .>= 0.0)
        @test all(mdrc.inc_p .>= 0.0)
        @test all(mdrc.inc_s .>= 0.0)

        # Check energy conservation
        P_energy, S_energy = energy_conservation(data)
        for i in eachindex(P_energy)
            @test isapprox(P_energy[i], 1.0, rtol=1e-2)
            @test isapprox(S_energy[i], 1.0, rtol=1e-2)
        end
        @info "P_energy: $(P_energy), S_energy: $(S_energy)"
    end

    @testset "Gaussian Surface - kspp ksmp support" begin
        # Create a Gaussian surface with small roughness
        σ = 10.0e-9  # RMS height (smaller for better numerical stability)
        ℓ = 100.0e-9  # Correlation length
        surface = GaussianSurface(σ, ℓ)

        # Create uniaxial medium with anisotropic properties
        ε_perp = -1.5 + 0.0im
        ε_para = -3.0 + 0.1im
        μ_perp = -0.5 + 0.0im
        μ_para = -0.5 + 0.0im
        below = Uniaxial(ε_perp, ε_para, μ_perp, μ_para)

        # More incident angles for a more comprehensive test
        paramconf = ParametersConfig(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)

        # Create solver data
        data = SolverData(paramconf, ensemble_iters, :reduced)
        params = data.params

        # Solve the ensemble problem
        solve_ensemble!(data)

        # Calculate MDRC
        mdrc = calc_mdrc(data)

        # Check that all reflectances and transmittances are non-negative
        @test all(mdrc.coh_p .>= 0.0)
        @test all(mdrc.coh_s .>= 0.0)
        @test all(mdrc.inc_p .>= 0.0)
        @test all(mdrc.inc_s .>= 0.0)

        # Check energy conservation
        P_energy, S_energy = energy_conservation(data)
        for i in eachindex(P_energy)
            @test isapprox(P_energy[i], 1.0, rtol=1e-2)
            @test isapprox(S_energy[i], 1.0, rtol=1e-2)
        end
        @info "P_energy: $(P_energy), S_energy: $(S_energy)"
    end
end
