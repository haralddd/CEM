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
sim_name = ARGS[1]
out_dir = "plots/$(sim_name)"
mkpath(out_dir)

metallic_dir = "output/$(sim_name)/metallic"
hyperbolic1_dir = "output/$(sim_name)/hyperbolic_type1"
hyperbolic2_dir = "output/$(sim_name)/hyperbolic_type2"

ylabel = L"\langle\partial R_{p}/\partial\theta_{s}\rangle_{\mathrm{incoh}}"
fig = Figure(fontsize=24, size=(600, 400))
ax = Axis(fig[1, 1],
    xlabel=L"\theta_s\ [\text{deg}]",
    ylabel=ylabel,
    xticks= -90:30:90,
    xminorticks = -90:10:90,
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorgridstyle = (:dot, :dense),
    xminorgridcolor = (:black, 0.12),
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
save(joinpath(out_dir, "metallic.pdf"), fig)

fig = Figure(fontsize=24, size=(600, 400))
ax = Axis(fig[1, 1],
    xlabel=L"\theta_s\ [\text{deg}]",
    ylabel=ylabel,
    xticks= -90:30:90,
    xminorticks = -90:10:90,
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorgridstyle = (:dot, :dense),
    xminorgridcolor = (:black, 0.12),
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
save(joinpath(out_dir, "hyperbolic_type1.pdf"), fig)

fig = Figure(fontsize=24, size=(600, 400))
ax = Axis(fig[1, 1],
    xlabel=L"\theta_s\ [\text{deg}]",
    ylabel=ylabel,
    xticks= -90:30:90,
    xminorticks = -90:10:90,
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorgridstyle = (:dot, :dense),
    xminorgridcolor = (:black, 0.12),
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
save(joinpath(out_dir, "hyperbolic_type2.pdf"), fig)