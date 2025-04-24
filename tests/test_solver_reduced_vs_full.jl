include("testconfig.jl")
using LinearAlgebra
using Statistics
using Test
using Random
using CairoMakie

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


    # Compare reduced solver with full solver for the same material
    @testset "Reduced vs Full Solver Comparison" begin
        # Create materials that both solvers can handle (μ_perp = μ_para)
        # Test with different permittivity values
        test_materials = [
            (2.0 + 0.0im, 2.0 + 0.0im),    # Isotropic case (εₚₑᵣₚ = εₚₐᵣₐ)
            (2.0 + 0.0im, 3.0 + 0.0im),    # Anisotropic case 1
            (-5.0 + 0.1im, -10.0 + 0.2im), # Anisotropic case 2 with complex values
            (1.5 + 0.0im, 4.0 + 0.0im)     # Anisotropic case 3
        ]

        # Test with Gaussian surface
        σ = 10.0e-9  # RMS height
        ℓ = 100.0e-9  # Correlation length
        surface = GaussianSurface(σ, ℓ)

        # Test with one representative material
        eps_perp = -17.0 + 0.48im
        eps_para = -17.0 + 0.48im
        below = Uniaxial(eps_perp, eps_para, μ_perp, μ_para)

        # Create parameters with multiple incident angles
        θs = [10.0]
        seed = 42
        params = Parameters(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs, seed=seed)

        # Create solver data for both solvers
        iters = 3  # Small number for test purposes
        data_reduced = SolverData(params, iters, :reduced)
        data_full = SolverData(params, iters, :full)

        # Set the same random seed for both solvers
        Random.seed!(seed)
        solve_ensemble!(data_reduced)

        # Reset the seed to ensure the same surfaces
        Random.seed!(seed)
        solve_ensemble!(data_full)

        # Calculate MDRC for both solvers
        mdrc_reduced = calc_mdrc(data_reduced)
        mdrc_full = calc_mdrc(data_full)

        err_coh_p = abs.(mdrc_reduced.coh_p - mdrc_full.coh_p)
        err_coh_s = abs.(mdrc_reduced.coh_s - mdrc_full.coh_s)
        err_inc_p = abs.(mdrc_reduced.inc_p - mdrc_full.inc_p)
        err_inc_s = abs.(mdrc_reduced.inc_s - mdrc_full.inc_s)


        @test all(isapprox.(err_coh_p, 0.0, atol=1e-15))
        @test all(isapprox.(err_coh_s, 0.0, atol=1e-15))
        @test all(isapprox.(err_inc_p, 0.0, atol=1e-15))
        @test all(isapprox.(err_inc_s, 0.0, atol=1e-15))

        fig = Figure()
        ax = fig[1, 1] = Axis(fig)
        lines!(ax, mdrc_full.θs, mdrc_full.inc_p[:], label="incoherent mdrc P full")
        lines!(ax, mdrc_full.θs, mdrc_full.inc_s[:], label="incoherent mdrc S full")
        lines!(ax, mdrc_reduced.θs, mdrc_reduced.inc_p[:], label="incoherent mdrc P reduced")
        lines!(ax, mdrc_reduced.θs, mdrc_reduced.inc_s[:], label="incoherent mdrc S reduced")
        axislegend(ax)
        display(fig)
        mkpath("plots/tests")
        save("plots/tests/test_reduced_vs_full.pdf", fig)
    end
end