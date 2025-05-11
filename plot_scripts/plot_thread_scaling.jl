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
fig = Figure(size=(600, 600))
ax1 = Axis(
    fig[1, 1],
    ylabel="Time (seconds)",
    title="Thread Scaling Performance",
    yscale=log10,
    xscale=log2
)
hidexdecorations!(ax1)

# Plot the data
reduced_line = lines!(ax1, threads, reduced_solver_times,
    color=:blue, linewidth=3,
    label="RRE")
scatter!(ax1, threads, reduced_solver_times,
    color=:blue, markersize=10)

full_line = lines!(ax1, threads, full_solver_times,
    color=:red, linewidth=3,
    label="Full")
scatter!(ax1, threads, full_solver_times,
    color=:red, markersize=10)

# Add thread count labels on x-axis
ax1.xticks = (threads, string.(threads))

# Calculate ideal scaling (T₁/n)
ideal_scaling_reduced = [reduced_solver_times[1] / n for n in threads]
ideal_scaling_full = [full_solver_times[1] / n for n in threads]

# Plot ideal scaling lines (dashed)
ideal_reduced = lines!(ax1, threads, ideal_scaling_reduced,
    color=:blue, linewidth=2, linestyle=:dash,
    label="Ideal (RRE)")

ideal_full = lines!(ax1, threads, ideal_scaling_full,
    color=:red, linewidth=2, linestyle=:dash,
    label="Ideal (Full)")

# Add grid lines
ax1.xgridvisible = true
ax1.ygridvisible = true

# Add legend
axislegend(ax1)

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
linkxaxes!(ax1, ax2)

# Plot speedup curves
lines!(ax2, threads, speedup_reduced,
    color=:blue, linewidth=3,
    label="RRE")
scatter!(ax2, threads, speedup_reduced,
    color=:blue, markersize=10)

lines!(ax2, threads, speedup_full,
    color=:red, linewidth=3,
    label="Full")
scatter!(ax2, threads, speedup_full,
    color=:red, markersize=10)

# Add ideal speedup line (dashed)
ideal_speedup = [n for n in threads]
lines!(ax2, threads, ideal_speedup,
    color=:black, linewidth=2, linestyle=:dash,
    label="Ideal")

# Add thread count labels on x-axis
ax2.xticks = (threads, string.(threads))

# Add grid lines
ax2.xgridvisible = true
ax2.ygridvisible = true

# Add legend
axislegend(ax2)

# Save the figure
save("thread_scaling.pdf", fig)