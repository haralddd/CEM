using RayleighSolver
using CairoMakie
using LaTeXStrings
using Colors

alpha = 1.0
linewidth = 2.0
markersize = 12
colors = [(:red, alpha), (:blue, alpha), (:green, alpha), (:orange, alpha)]

metallic_dir = "output/sweep/metallic"

fig = Figure(fontsize=24, size=(800, 600))
ax = Axis(fig[1, 1],
    xlabel=L"\theta_s\ [^\circ]",
    ylabel=L"\text{Inc. MDRC}",
    xticks= -90:30:90,
    xgridvisible=false,
    ygridvisible=false
)

for (i, fname) in enumerate(readdir(metallic_dir))
    label = splitext(fname)[1]
    data = load_solver_data(joinpath(metallic_dir, fname))
    mdrc_data = calc_mdrc(data)
    
    ln = lines!(ax, 
        mdrc_data.θs, mdrc_data.inc_p[:, 1], 
        label=L"\varepsilon_\perp = %$(label) \varepsilon",
        color=colors[i],
        linewidth=linewidth)

    ymax, maxidx = findmax(mdrc_data.inc_p[:, 1])
    sc = scatter!(ax, mdrc_data.θs[maxidx], ymax, color=colors[i], marker=:hline, markersize=markersize)

    translate!(sc, 0, 0, 10000-ymax)
    translate!(ln, 0, 0, 10000-ymax)
end
axislegend(ax)
save("plots/metallic.pdf", fig)

hyperbolic1_dir = "output/sweep/hyperbolic_type1"

fig = Figure(fontsize=24, size=(800, 600))
ax = Axis(fig[1, 1],
    xlabel=L"\theta_s\ [^\circ]",
    ylabel=L"\text{Inc. MDRC}",
    xticks= -90:30:90,
    xgridvisible=false,
    ygridvisible=false
)

for (i, fname) in enumerate(readdir(hyperbolic1_dir))
    label = splitext(fname)[1]
    data = load_solver_data(joinpath(hyperbolic1_dir, fname))
    mdrc_data = calc_mdrc(data)
    
    ln = lines!(ax, 
        mdrc_data.θs, mdrc_data.inc_p[:, 1], 
        label=L"\varepsilon_\parallel = %$(label) \varepsilon",
        color=colors[i],
        linewidth=linewidth)
    ymax, maxidx = findmax(mdrc_data.inc_p[:, 1])
    translate!(ln, 0, 0, 10000-ymax)
end
axislegend(ax)
save("plots/hyberbolic_type1.pdf", fig)

hyperbolic2_dir = "output/sweep/hyperbolic_type2"

fig = Figure(fontsize=24, size=(800, 600))
ax = Axis(fig[1, 1],
    xlabel=L"\theta_s\ [^\circ]",
    ylabel=L"\text{Inc. MDRC}",
    xticks= -90:30:90,
    xgridvisible=false,
    ygridvisible=false
)

for (i, fname) in enumerate(readdir(hyperbolic2_dir))
    label = splitext(fname)[1]
    data = load_solver_data(joinpath(hyperbolic2_dir, fname))
    mdrc_data = calc_mdrc(data)
    
    ln = lines!(ax, 
        mdrc_data.θs, mdrc_data.inc_p[:, 1], 
        label=L"\varepsilon_\perp = %$(label) \varepsilon",
        color=colors[i],
        linewidth=linewidth)

    ymax, maxidx = findmax(mdrc_data.inc_p[:, 1])
    sc = scatter!(ax, mdrc_data.θs[maxidx], ymax, color=colors[i], marker=:hline, markersize=markersize)
    translate!(sc, 0, 0, 10000-ymax)
    translate!(ln, 0, 0, 10000-ymax)
end
axislegend(ax)
save("plots/hyberbolic_type2.pdf", fig)