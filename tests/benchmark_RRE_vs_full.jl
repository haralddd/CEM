using RayleighSolver
using BenchmarkTools
using Random
using Statistics
using Printf
using CairoMakie

include("test_RT_relation.jl")  # Include the file with calc_Tp_Ts function

"""
    benchmark_cem_methods(Nx_values, repeats=3)

Compare timing between solving the Reduced Rayleigh Equations (using calc_Tp_Ts) 
and solving the full scattering system for different system sizes.

# Arguments
- `Nx_values`: Array of system sizes to test
- `repeats`: Number of repetitions for each benchmark to ensure accurate timing

# Returns
- `results`: Dictionary with timing results
"""
function benchmark_cem_methods(Nx_values, repeats=3)
    results = Dict(
        :Nx => Nx_values,
        :rre_times => Float64[],
        :full_times => Float64[]
    )
    
    # Using a simple system configuration for benchmarking
    θ0 = [0.0]
    Lx = 200.0
    tio2 = Uniaxial(8.43 + 0.0im, 6.84 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im)
    
    for Nx in Nx_values
        println("\n=== Testing system size Nx = $Nx ===")
        
        # Setup parameters for reduced solver
        params_reduced = ParametersConfig(
            lambda = 589.3e-9,
            Lx = Lx,
            Nx = Nx,
            θs = θ0,
            below = tio2,
            solver_type = :reduced
        )
        
        # Setup parameters for full solver
        params_full = ParametersConfig(
            lambda = 589.3e-9,
            Lx = Lx,
            Nx = Nx,
            θs = θ0,
            below = tio2,
            solver_type = :full
        )
        
        # Generate surface (same for both methods)
        surf = GaussianSurface(10.0e-9, 100.0e-9)
        
        # Create solver data objects
        data_reduced = SolverData(params_reduced)
        data_full = SolverData(params_full)
        
        # Create allocation and precompute objects
        pre_reduced = Precomputed(data_reduced)
        alloc_reduced = Preallocated(data_reduced)
        
        pre_full = Precomputed(data_full)
        alloc_full = Preallocated(data_full)
        
        # Generate same surface for both methods
        Random.seed!(12345)
        generate_surface!(alloc_reduced, surf)
        Random.seed!(12345)
        generate_surface!(alloc_full, surf)
        
        # Precompute values
        precompute!(pre_reduced, data_reduced)
        precompute!(pre_full, data_full)
        
        # Benchmark Reduced Rayleigh Equations (RRE) method with calc_Tp_Ts
        rre_times = Float64[]
        for _ in 1:repeats
            # Solve RRE to get R
            t_start = time()
            solve_single!(alloc_reduced, pre_reduced, data_reduced)
            
            # Use R to calculate T via calc_Tp_Ts
            Tp = zeros(ComplexF64, size(alloc_reduced.Rp))
            Ts = zeros(ComplexF64, size(alloc_reduced.Rs))
            calc_Tp_Ts!(Tp, Ts, alloc_reduced.Rp, alloc_reduced.Rs, pre_reduced, data_reduced.params)
            t_end = time()
            
            push!(rre_times, t_end - t_start)
        end
        
        # Benchmark Full Solver method
        full_times = Float64[]
        for _ in 1:repeats
            t_start = time()
            solve_single!(alloc_full, pre_full, data_full)
            t_end = time()
            push!(full_times, t_end - t_start)
        end
        
        # Store median times
        push!(results[:rre_times], median(rre_times))
        push!(results[:full_times], median(full_times))
        
        # Print results
        @printf("RRE method:  %.4f seconds\n", median(rre_times))
        @printf("Full method: %.4f seconds\n", median(full_times))
        @printf("Speedup: %.2fx\n", median(full_times) / median(rre_times))
    end
    
    return results
end

"""
    plot_benchmark_results(results)

Plot the benchmark results comparing RRE vs Full solver performance.
"""
function plot_benchmark_results(results)
    Nx_values = results[:Nx]
    rre_times = results[:rre_times]
    full_times = results[:full_times]
    speedups = full_times ./ rre_times
    
    # Create figure with two subplots
    fig = Figure(resolution=(1000, 800))
    
    # Plot execution times
    ax1 = Axis(fig[1, 1], 
        xlabel = "System size (Nx)",
        ylabel = "Execution time (seconds)",
        title = "RRE vs Full solver performance",
        yscale = log10,
        xscale = log2
    )
    
    # Add reference lines for O(N log N) and O(N^2) complexity
    Nx_range = range(minimum(Nx_values), maximum(Nx_values), length=100)
    n_log_n_scale = rre_times[1] * (Nx_range .* log2.(Nx_range)) ./ (Nx_values[1] * log2(Nx_values[1]))
    n_squared_scale = full_times[1] * (Nx_range.^2) ./ (Nx_values[1]^2)
    
    lines!(ax1, Nx_range, n_log_n_scale, linestyle=:dash, color=:blue, linewidth=1.5, 
           label="O(N log N)")
    lines!(ax1, Nx_range, n_squared_scale, linestyle=:dash, color=:red, linewidth=1.5, 
           label="O(N²)")
    
    scatterlines!(ax1, Nx_values, rre_times, 
        label = "RRE + calc_Tp_Ts", 
        color = :blue, 
        marker = :circle,
        markersize = 10)
        
    scatterlines!(ax1, Nx_values, full_times, 
        label = "Full solver", 
        color = :red, 
        marker = :square,
        markersize = 10)
    
    ax1.xticks = (Nx_values, string.(Nx_values))
    axislegend(ax1, position=:lt)
    
    # Plot speedups
    ax2 = Axis(fig[2, 1], 
        xlabel = "System size (Nx)",
        ylabel = "Speedup factor (Full / RRE)",
        title = "Speedup of RRE method vs Full solver",
        xscale = log2
    )
    
    scatterlines!(ax2, Nx_values, speedups, 
        color = :green, 
        marker = :diamond,
        markersize = 10)
    
    # Add reference line for O(N) speedup
    n_scale = speedups[1] * (Nx_range) ./ (Nx_values[1])
    lines!(ax2, Nx_range, n_scale, linestyle=:dash, color=:green, linewidth=1.5, 
           label="O(N)")
           
    ax2.xticks = (Nx_values, string.(Nx_values))
    axislegend(ax2)
        
    # Save the figure
    save("benchmark_rre_vs_full.pdf", fig)
    
    return fig
end

# Test different system sizes (using powers of 2 for better FFT performance)
Nx_values = [512, 1024, 2048, 4096, 8192]
results = benchmark_cem_methods(Nx_values)

# Plot results
fig = plot_benchmark_results(results)
display(fig)

# Print table of results
println("\n=== Summary of Results ===")
println("System Size | RRE Time (s) | Full Time (s) | Speedup")
println("-----------|-------------|-------------|--------")
for i in eachindex(results[:Nx])
    @printf("%11d | %11.4f | %11.4f | %7.2fx\n", 
            results[:Nx][i], 
            results[:rre_times][i], 
            results[:full_times][i], 
            results[:full_times][i] / results[:rre_times][i])
end
