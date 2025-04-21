include("testconfig.jl")
include("test_surface.jl")
using Statistics
using Peaks
using CairoMakie
using LaTeXStrings
using DataFrames

# Enable LaTeX rendering in CairoMakie
CairoMakie.activate!(type = "pdf")

function calculate_errors(data::SolverData, iters::Int)
    params = data.params
    alloc = Preallocated(data)
    surf = params.surf
    xs = params.xs
    dx = params.dx
    d = params.surf.d
    
    # Analytical values
    analytical_mean = 0.0
    analytical_s = analytical_slope(surf)
    analytical_D = analytical_dist(surf)
    analytical_rms = d
    
    # Computed values
    meanval = 0.0
    rms = 0.0
    slope = 0.0
    dist = 0.0
    
    for m in 1:iters
        generate_surface!(alloc, params)
        meanval = observe(meanval, mean(alloc.ys), m)
        rms = observe(rms, sqrt(mean(alloc.ys .^ 2)), m)
        slope = observe(slope, mean_slope(alloc.ys, dx), m)
        dist = observe(dist, mean_peak_valley_dist(xs, alloc.ys), m)
    end
    
    # Calculate relative errors
    mean_error = abs(meanval - analytical_mean)
    rms_error = abs(rms - analytical_rms) / analytical_rms
    slope_error = abs(slope - analytical_s) / analytical_s
    dist_error = abs(dist - analytical_D) / analytical_D
    
    return Dict(
        "dx" => dx,
        "Nx" => params.Nx,
        "mean_error" => mean_error,
        "rms_error" => rms_error,
        "slope_error" => slope_error,
        "dist_error" => dist_error
    )
end

function test_resolution_scaling(surf_type::Symbol, Nx_values, Lx, lambda, iters=1000)
    results = DataFrame(
        dx = Float64[],
        Nx = Int[],
        mean_error = Float64[],
        rms_error = Float64[],
        slope_error = Float64[],
        dist_error = Float64[]
    )

    σ = 10.0e-9
    
    for Nx in Nx_values
        params = if surf_type == :gaussian
            Parameters(Nx=Nx, Lx=Lx, lambda=lambda, surf=GaussianSurface(σ, 100.0e-9))
        elseif surf_type == :rectangular
            Parameters(Nx=Nx, Lx=Lx, lambda=lambda, surf=RectangularSurface(σ, 0.70, 1.30))
        else
            error("Unknown surface type: $surf_type")
        end
        
        data = SolverData(params, iters)
        errors = calculate_errors(data, iters)
        push!(results, errors)
        
        # Print progress
        println("Completed Nx = $Nx for $surf_type surface")
    end
    
    return results
end

function plot_scaling_results(gaussian_results, rect_results)
    # Separate plots for each error metric
    plot_metrics = ["mean_error", "rms_error", "slope_error", "dist_error"]
    plot_titles = [L"\mathrm{Error\ in\ mean\ value}", L"\mathrm{Error\ in\ RMS\ height}", L"\mathrm{Error\ in\ mean\ slope}", L"\mathrm{Error\ in\ peak-valley\ distance}"]

    # Save plots
    mkpath("plots/resolution_scaling")
    for (metric, title) in zip(plot_metrics, plot_titles)
        fig = Figure(size = (800, 500), fontsize=24)
        ax = Axis(fig[1, 1],
            title = title,
            xlabel = L"N_x",
            ylabel = L"$\mathrm{Relative\ Error} \ [%]$",
            xticks = [2048, 4096, 6144, 8192])
        
        # For Gaussian surface
        scatter!(ax, gaussian_results.Nx, gaussian_results[:, metric] * 100, 
                 label = L"\mathrm{Gaussian}", marker = :circle, markersize = 12, color = :steelblue)
        
        # For Rectangular surface
        scatter!(ax, rect_results.Nx, rect_results[:, metric] * 100, 
                 label = L"\mathrm{West-O'Donnell}", marker = :rect, markersize = 12, color = :crimson)
        
        ylims!(ax, 0, nothing)
        # Add legend
        axislegend(ax, position = :rt)
        save("plots/resolution_scaling/$(metric)_scaling.pdf", fig)

    end
    
    return nothing
end

# Run the scaling tests
function run_scaling_tests()
    # Create plots directory if it doesn't exist
    mkpath("plots/resolution_scaling")
    
    # Define parameters for the tests
    lambda = 632.8e-9
    Lx = 100 * lambda
    k0 = 2π / lambda
    
    Nx_values = range(2048, 8192, step=1024)
    @info "Q low: $(π / (Lx / Nx_values[1]) / k0)"
    @info "Q high: $(π / (Lx / Nx_values[end]) / k0)"
    
    # Lower iterations for faster testing (adjust based on your needs)
    iters = 2000
    
    println("Running Gaussian surface resolution tests...")
    gaussian_results = test_resolution_scaling(:gaussian, Nx_values, Lx, lambda, iters)
    
    println("Running Rectangular surface resolution tests...")
    rect_results = test_resolution_scaling(:rectangular, Nx_values, Lx, lambda, iters)
    
    # Plot the results
    println("Plotting results...")
    plot_scaling_results(gaussian_results, rect_results)
    
    println("Results saved to plots/ directory")
    return gaussian_results, rect_results
end

# Execute tests
gaussian_results, rect_results = run_scaling_tests()

# Print the data tables
println("\nGaussian Surface Results:")
println(gaussian_results)

println("\nRectangular Surface Results:")
println(rect_results)
