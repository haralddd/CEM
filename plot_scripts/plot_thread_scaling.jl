using CairoMakie
using LaTeXStrings

# Thread counts
threads = [1, 2, 4, 8, 16, 32]
reduced_solver_times = [
    # Replace these values with your actual times for reduced solver
    1857.204560, 938.639220, 475.763184, 256.642851, 147.645244, 144.869521
]

full_solver_times = [
    # Replace these values with your actual times for full solver
    13105.974891, 6549.320591, 3310.622496, 1686.997213, 921.704879, 912.475211
]

# Create the plot
fig = Figure(resolution=(600, 400))
ax = Axis(
    fig[1, 1],
    xlabel="Number of Threads",
    ylabel="Time (seconds)",
    title="Thread Scaling Performance",
    yscale=log10,
    xscale=log2
)

# Plot the data
reduced_line = lines!(ax, threads, reduced_solver_times, 
                     color=:blue, linewidth=3, 
                     label="Reduced Solver")
scatter!(ax, threads, reduced_solver_times, 
        color=:blue, markersize=10)

full_line = lines!(ax, threads, full_solver_times, 
                  color=:red, linewidth=3, 
                  label="Full Solver")
scatter!(ax, threads, full_solver_times, 
        color=:red, markersize=10)

# Add thread count labels on x-axis
ax.xticks = (threads, string.(threads))

# Calculate ideal scaling (T₁/n)
ideal_scaling_reduced = [reduced_solver_times[1] / n for n in threads]
ideal_scaling_full = [full_solver_times[1] / n for n in threads]

# Plot ideal scaling lines (dashed)
ideal_reduced = lines!(ax, threads, ideal_scaling_reduced, 
                     color=:blue, linewidth=2, linestyle=:dash, 
                     label="Ideal Scaling (Reduced)")

ideal_full = lines!(ax, threads, ideal_scaling_full, 
                  color=:red, linewidth=2, linestyle=:dash, 
                  label="Ideal Scaling (Full)")

# Add grid lines
ax.xgridvisible = true
ax.ygridvisible = true

# Add legend
Legend(fig[1, 2], [reduced_line, full_line, ideal_reduced, ideal_full], 
      ["Reduced Solver", "Full Solver", "Ideal Scaling (Reduced)", "Ideal Scaling (Full)"])

# Calculate speedup
speedup_reduced = [reduced_solver_times[1] / t for t in reduced_solver_times]
speedup_full = [full_solver_times[1] / t for t in full_solver_times]

# Create a second y-axis for speedup values
ax2 = Axis(
    fig[2, 1],
    xlabel="Number of Threads",
    ylabel="Speedup (T₁/Tₙ)",
    xscale=log2
)

# Plot speedup curves
lines!(ax2, threads, speedup_reduced, 
      color=:blue, linewidth=3, 
      label="Reduced Solver")
scatter!(ax2, threads, speedup_reduced, 
        color=:blue, markersize=10)

lines!(ax2, threads, speedup_full, 
      color=:red, linewidth=3, 
      label="Full Solver")
scatter!(ax2, threads, speedup_full, 
        color=:red, markersize=10)

# Add ideal speedup line (dashed)
ideal_speedup = [n for n in threads]
lines!(ax2, threads, ideal_speedup, 
      color=:black, linewidth=2, linestyle=:dash, 
      label="Ideal Speedup")

# Add thread count labels on x-axis
ax2.xticks = (threads, string.(threads))

# Add grid lines
ax2.xgridvisible = true
ax2.ygridvisible = true

# Add legend
Legend(fig[2, 2], 
      ["Reduced Solver", "Full Solver", "Ideal Speedup"], 
      ["Reduced Solver", "Full Solver", "Ideal Speedup"])

# Save the figure
save("thread_scaling.pdf", fig)
println("Plot saved as 'thread_scaling.pdf'")

# Print a table of results
println("\nThread Scaling Results:")
println("Threads | Reduced Time (s) | Reduced Speedup | Full Time (s) | Full Speedup")
println("--------|-----------------|----------------|--------------|-------------")
for i in 1:length(threads)
    println("$(threads[i]) | $(round(reduced_solver_times[i], digits=2)) | $(round(speedup_reduced[i], digits=2)) | $(round(full_solver_times[i], digits=2)) | $(round(speedup_full[i], digits=2))")
end

# Create a parser function for extracting benchmark data from output files
"""
    parse_benchmark_output(filename::String)

Parse the output file from thread scaling benchmarks.
Returns the times for reduced and full solvers as vectors.
"""
function parse_benchmark_output(filename::String)
    reduced_times = Float64[]
    full_times = Float64[]
    
    current_section = :none
    
    open(filename) do f
        for line in eachline(f)
            if occursin("Default uniaxial reduced ensemble solve", line)
                current_section = :reduced
                continue
            elseif occursin("Default uniaxial full ensemble solve", line)
                current_section = :full
                continue
            end
            
            if current_section != :none && occursin("seconds", line)
                # Extract timing information
                time_match = match(r"([0-9]+\.[0-9]+) seconds", line)
                if time_match !== nothing
                    time_value = parse(Float64, time_match.captures[1])
                    
                    if current_section == :reduced
                        push!(reduced_times, time_value)
                    elseif current_section == :full
                        push!(full_times, time_value)
                    end
                end
            end
        end
    end
    
    return reduced_times, full_times
end

# Example usage:
# reduced_times, full_times = parse_benchmark_output("path/to/output_file.out")
# Then replace the sample data with these parsed times
