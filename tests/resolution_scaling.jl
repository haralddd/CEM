include("testconfig.jl")
include("test_surface.jl")
using Statistics
using Peaks
using CairoMakie
using LaTeXStrings
using DataFrames
using FFTW

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
    rmsval = 0.0
    slope = 0.0
    dist = 0.0
    
    for m in 1:iters
        generate_surface!(alloc, params)
        cur_mean = mean(alloc.ys)
        meanval = observe(meanval, cur_mean, m)
        rmsval = observe(rmsval, rms(alloc.ys, cur_mean), m)
        slope = observe(slope, rms_slope(alloc.ys, dx), m)
        dist = observe(dist, mean_peak_valley_dist(xs, alloc.ys), m)
    end
    
    # Calculate relative errors
    mean_error = abs(meanval - analytical_mean)
    rms_error = abs(rmsval - analytical_rms) / analytical_rms
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
            ParametersConfig(Nx=Nx, Lx=Lx, lambda=lambda, surf=GaussianSurface(σ, 100.0e-9))
        elseif surf_type == :rectangular
            ParametersConfig(Nx=Nx, Lx=Lx, lambda=lambda, surf=RectangularSurface(σ, 0.70, 1.30))
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
    plot_titles = [L"\mathrm{Error\ in\ mean\ value}", L"\mathrm{Error\ in\ RMS\ height}", L"\mathrm{Error\ in\ RMS\ slope}", L"\mathrm{Error\ in\ mean\ peak-valley\ distance}"]

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

# Calculate power spectrum from surface realization
function calculate_power_spectrum(data::SolverData)
    # Get the surface data and parameters
    params = data.params
    alloc = Preallocated(data)
    generate_surface!(alloc, params)
    
    ys = copy(alloc.ys)
    dx = params.dx
    Nx = params.Nx
    
    # Remove mean to avoid DC component dominance
    ys .-= mean(ys)
    
    # Calculate FFT
    fft_result = fft(ys)
    
    # Calculate wavenumbers
    dk = 2π / (Nx * dx)  # Wavenumber resolution
    kx = fftshift(fftfreq(Nx, 1.0/dx))  # Shifted wavenumbers
    
    # Calculate power spectrum (normalize for numerical comparison)
    power_spectrum = abs2.(fftshift(fft_result)) .* dx / (2π * params.surf.d^2)
    
    # Only return the positive wavenumbers (symmetry in power spectrum)
    pos_indices = kx .>= 0
    return kx[pos_indices], power_spectrum[pos_indices]
end

# Calculate analytical power spectrum based on surface type
function analytical_power_spectrum(kx::Vector{Float64}, surf::RandomSurface)
    # Calculate analytical power spectrum based on the correlation function
    # The correlation function is related to the power spectrum
    if isa(surf, GaussianSurface)
        # For Gaussian surfaces
        a = surf.a
        return [√π * a * exp(-(a * k / 2.0)^2) for k in kx]
    elseif isa(surf, RectangularSurface)
        # For Rectangular surfaces
        kp = surf.kp
        km = surf.km
        return [if (km <= abs(k) <= kp) π / (kp - km) else 0.0 end for k in kx]
    else
        error("Unknown surface type: $(typeof(surf))")
    end
end

# Plot power spectrum for a given surface
function plot_power_spectrum(surf_type::Symbol, Nx::Int, Lx::Float64, lambda::Float64, iters::Int)
    # Set up parameters
    σ = 10.0e-9  # RMS height
    
    params = if surf_type == :gaussian
        ParametersConfig(Nx=Nx, Lx=Lx, lambda=lambda, surf=GaussianSurface(σ, 100.0e-9))
    elseif surf_type == :rectangular
        ParametersConfig(Nx=Nx, Lx=Lx, lambda=lambda, surf=RectangularSurface(σ, 0.70, 1.30))
    else
        error("Unknown surface type: $surf_type")
    end
    
    data = SolverData(params, iters)
    
    # Calculate wavenumbers and initialize average power spectrum
    kx, power = calculate_power_spectrum(data)
    avg_power = zeros(length(power))
    
    # Average over multiple realizations
    for i in 1:iters
        _, power = calculate_power_spectrum(data)
        avg_power .+= power
    end
    avg_power ./= iters
    
    # Calculate analytical power spectrum
    analytical = analytical_power_spectrum(kx, data.params.surf)
    
    # Create plot
    fig = Figure(size = (800, 500), fontsize=24)
    ax = Axis(fig[1, 1],
        title = L"\mathrm{Power\ Spectrum}",
        xlabel = L"k\ [\mathrm{rad/m}]",
        ylabel = L"g(k)",
        xscale = log10,
        yscale = log10)
    
    # Normalize to k0 for better readability
    k0 = 2π / lambda
    
    # Plot numerical results
    lines!(ax, kx/k0, avg_power, linewidth=2, label=L"\mathrm{Numerical}", color=:steelblue)
    
    # Plot analytical results
    lines!(ax, kx/k0, analytical, linewidth=2, label=L"\mathrm{Analytical}", color=:crimson, linestyle=:dash)
    
    # Add legend
    axislegend(ax, position = :rt)
    
    # Determine surface type label for filename
    surf_label = surf_type == :gaussian ? "gaussian" : "rectangular"
    
    # Save plot
    save_path = "plots/power_spectrum_$(surf_label)_Nx$(Nx).pdf"
    save(save_path, fig)
    
    return kx, avg_power, analytical
end

# Run the scaling tests
function run_scaling_tests()
    # Create plots directory if it doesn't exist
    mkpath("plots/resolution_scaling")
    mkpath("plots/power_spectrum")
    
    # Define parameters for the tests
    lambda = 632.8e-9
    Lx = 100
    k0 = 2π / lambda
    
    Nx_values = range(2048, 8192, step=1024)
    @info "Q low: $(π / (Lx * lambda / Nx_values[1]) / k0)"
    @info "Q high: $(π / (Lx * lambda / Nx_values[end]) / k0)"
    
    # Lower iterations for faster testing (adjust based on your needs)
    iters = 2000
    power_iters = 100  # Use fewer iterations for power spectrum calculations
    
    println("Running Gaussian surface resolution tests...")
    gaussian_results = test_resolution_scaling(:gaussian, Nx_values, Lx, lambda, iters)
    
    println("Running Rectangular surface resolution tests...")
    rect_results = test_resolution_scaling(:rectangular, Nx_values, Lx, lambda, iters)
    
    # Plot the results
    println("Plotting resolution scaling results...")
    plot_scaling_results(gaussian_results, rect_results)
    
    # Calculate and plot power spectra for highest resolution
    println("\nCalculating and plotting power spectra...")
    highest_Nx = Int(Nx_values[end])
    
    println("  - Gaussian surface power spectrum...")
    gaussian_kx, gaussian_power, gaussian_analytical = plot_power_spectrum(:gaussian, highest_Nx, Float64(Lx), lambda, power_iters)
    
    println("  - Rectangular surface power spectrum...")
    rect_kx, rect_power, rect_analytical = plot_power_spectrum(:rectangular, highest_Nx, Float64(Lx), lambda, power_iters)
    
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
