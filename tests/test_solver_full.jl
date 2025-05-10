using RayleighSolver
using Test

@testset "Full Uniaxial Solver Tests" begin
    # Common parameters
    θ0 = 20.0  # default incident angle in degrees
    λ = 632.8  # wavelength in nm
    
    # Surface parameters
    Nx = 4096  # Number of discretization points
    Lx = 100  # Length of the surface
    
    # Media
    above = Vacuum()
    ensemble_iters = 1000
    θs = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]
    mu = 1.0 + 0.0im
    
    @testset "Flat Surface - Transparent Anisotropic Permeability" begin
        # Create uniaxial medium with anisotropic permeability (μₚₑᵣₚ ≠ μₚₐᵣₐ)
        # This is valid for the full solver but not for the reduced solver
        ε_perp = 2.0 + 0.0im
        ε_para = 2.0 + 0.0im
        
        # Create flat surface and uniaxial medium
        surface = FlatSurface()
        below = Uniaxial(ε_perp, ε_para, mu, mu)
        
        # Create parameters and solver data
        paramconf = ParametersConfig(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)
        data = SolverData(paramconf, 1, :full)
        params = data.params
        
        # Preallocate and precompute
        precomputed = Precomputed(data)
        precompute!(precomputed, data)
        
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
    
    @testset "Gaussian Surface - Bianisotropy Transparent" begin
        # Create a Gaussian surface with small roughness
        σ = 10.0e-9  # RMS height (smaller for better numerical stability)
        ℓ = 100.0e-9  # Correlation length
        surface = GaussianSurface(σ, ℓ)
        
        # Create uniaxial medium with anisotropic properties
        ε_perp = 2.0 + 0.0im
        ε_para = 3.0 + 0.0im
        below = Uniaxial(ε_perp, ε_para, mu, mu)
        
        # More incident angles for a more comprehensive test
        paramconf = ParametersConfig(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)
        
        # Create solver data
        data = SolverData(paramconf, ensemble_iters, :full)
        params = data.params
        
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

    @testset "Gaussian Surface - Bianisotropy HMM1" begin
        # Create a Gaussian surface with small roughness
        σ = 10.0e-9  # RMS height (smaller for better numerical stability)
        ℓ = 100.0e-9  # Correlation length
        surface = GaussianSurface(σ, ℓ)

        # Create uniaxial medium with anisotropic properties
        ε_perp = -7.0 + 0.0im
        ε_para = 3.0 + 0.0im
        below = Uniaxial(ε_perp, ε_para, mu, mu)

        # More incident angles for a more comprehensive test
        paramconf = ParametersConfig(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)

        # Create solver data
        data = SolverData(paramconf, ensemble_iters, :full)
        params = data.params

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

    @testset "Gaussian Surface - Bianisotropy HMM2" begin
        # Create a Gaussian surface with small roughness
        σ = 10.0e-9  # RMS height (smaller for better numerical stability)
        ℓ = 100.0e-9  # Correlation length
        surface = GaussianSurface(σ, ℓ)

        # Create uniaxial medium with anisotropic properties
        ε_perp = 3.0 + 0.0im
        ε_para = -7.0 + 0.0im
        below = Uniaxial(ε_perp, ε_para, mu, mu)

        # More incident angles for a more comprehensive test
        paramconf = ParametersConfig(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)

        # Create solver data
        data = SolverData(paramconf, ensemble_iters, :full)
        params = data.params

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
        ε_para = -3.0 + 0.01im
        μ_perp = -0.5 + 0.0im
        μ_para = -0.5 + 0.0im
        below = Uniaxial(ε_perp, ε_para, μ_perp, μ_para)

        # More incident angles for a more comprehensive test
        paramconf = ParametersConfig(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)

        # Create solver data
        data = SolverData(paramconf, ensemble_iters, :full)
        params = data.params

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

        # Check energy conservation
        P_energy, S_energy = energy_conservation(data)
        for i in eachindex(P_energy)
            @test isapprox(P_energy[i], 1.0, rtol=1e-2)
            @test isapprox(S_energy[i], 1.0, rtol=1e-2)
        end
        @info "P_energy: $(P_energy), S_energy: $(S_energy)"
    end
end
