using RayleighSolver
using CairoMakie
using LaTeXStrings
using Colors

alpha = 1.0
linewidth = 1.0
markersize = 10
colors = [:darkgreen, :royalblue, :firebrick, :darkorange]
markers = [:circle, :rect, :diamond, :utriangle]

axis_label_fontsize=32
metallic_dir = "output/sweep/metallic"

ylabel = L"\langle\partial R_{p}/\partial\theta_{s}\rangle_{\mathrm{incoh}}"
fig = Figure(fontsize=24, size=(600, 400))
ax = Axis(fig[1, 1],
    xlabel=L"\theta_s\ [\text{deg}]",
    ylabel=ylabel,
    xticks= -90:30:90,
    xgridvisible=false,
    ygridvisible=false,
    xlabelsize=axis_label_fontsize,
    ylabelsize=axis_label_fontsize
)

for (i, fname) in enumerate(readdir(metallic_dir))
    label = splitext(fname)[1]
    data = load_solver_data(joinpath(metallic_dir, fname))
    mdrc_data = calc_mdrc(data)
    
    ln = lines!(ax, 
        mdrc_data.θs, mdrc_data.inc_p[:, 1], 
        color=colors[i],
        linewidth=linewidth)

    ymax, maxidx = findmax(mdrc_data.inc_p[:, 1])
    sc = scatterlines!(ax, mdrc_data.θs[maxidx], ymax, color=colors[i], marker=markers[i], markersize=markersize,
        label=L"\varepsilon_\perp' = %$(label) \varepsilon'")

    translate!(sc, 0, 0, 10000-ymax)
    translate!(ln, 0, 0, 10000-ymax)
end
axislegend(ax)
save("plots/metallic.pdf", fig)

hyperbolic1_dir = "output/sweep/hyperbolic_type1"

fig = Figure(fontsize=24, size=(600, 400))
ax = Axis(fig[1, 1],
    xlabel=L"\theta_s\ [\text{deg}]",
    ylabel=ylabel,
    xticks= -90:30:90,
    xgridvisible=false,
    ygridvisible=false,
    xlabelsize=axis_label_fontsize,
    ylabelsize=axis_label_fontsize
)

for (i, fname) in enumerate(readdir(hyperbolic1_dir))
    label = splitext(fname)[1]
    if label == "1.0"
        continue
    end
    data = load_solver_data(joinpath(hyperbolic1_dir, fname))
    mdrc_data = calc_mdrc(data)
    
    ln = lines!(ax, 
        mdrc_data.θs, mdrc_data.inc_p[:, 1], 
        color=colors[i],
        linewidth=linewidth)
    ymax, maxidx = findmax(mdrc_data.inc_p[:, 1])
    sc = scatterlines!(ax, mdrc_data.θs[maxidx], ymax, color=colors[i], marker=markers[i], markersize=markersize,
        label=L"\varepsilon_\parallel' = %$(label) \varepsilon'")
    translate!(ln, 0, 0, 10000-ymax)
end
axislegend(ax)
save("plots/hyberbolic_type1.pdf", fig)

hyperbolic2_dir = "output/sweep/hyperbolic_type2"

fig = Figure(fontsize=24, size=(600, 400))
ax = Axis(fig[1, 1],
    xlabel=L"\theta_s\ [\text{deg}]",
    ylabel=ylabel,
    xticks= -90:30:90,
    xgridvisible=false,
    ygridvisible=false,
    xlabelsize=axis_label_fontsize,
    ylabelsize=axis_label_fontsize
)

for (i, fname) in enumerate(readdir(hyperbolic2_dir))
    label = splitext(fname)[1]
    if label == "1.0"
        continue
    end
    data = load_solver_data(joinpath(hyperbolic2_dir, fname))
    mdrc_data = calc_mdrc(data)
    
    ln = lines!(ax, 
        mdrc_data.θs, mdrc_data.inc_p[:, 1],
        color=colors[i],
        linewidth=linewidth)

    ymax, maxidx = findmax(mdrc_data.inc_p[:, 1])
    sc = scatterlines!(ax, mdrc_data.θs[maxidx], ymax, color=colors[i], marker=markers[i], markersize=markersize,
        label=L"\varepsilon_\perp' = %$(label) \varepsilon'")
    translate!(sc, 0, 0, 10000-ymax)
    translate!(ln, 0, 0, 10000-ymax)
end
axislegend(ax)
save("plots/hyberbolic_type2.pdf", fig)