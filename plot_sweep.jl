using RayleighSolver
using CairoMakie
using LaTeXStrings

metallic_dir = "output/sweep/metallic"

fig = Figure(fontsize=24, size=(800, 600))
ax = Axis(fig[1, 1],
    title = "MDRC",
    xlabel=L"\theta_s [deg]",
    ylabel=L"\text{MDRC}",
    xticks= -90:30:90,
    xgridvisible=false,
    ygridvisible=false
)

for fname in readdir(metallic_dir)
    label = splitext(fname)[1]
    data = load_solver_data(joinpath(metallic_dir, fname))
    mdrc_data = calc_mdrc(data)
    
    lines!(ax, mdrc_data.θs, mdrc_data.inc_p[:, 1], label=L"a = %$(label)")
end
axislegend(ax)
save("plots/metallic.pdf", fig)

hyperbolic1_dir = "output/sweep/hyperbolic_type1"

fig = Figure(fontsize=24, size=(800, 600))
ax = Axis(fig[1, 1],
    title = "MDRC",
    xlabel=L"\theta_s [deg]",
    ylabel=L"\text{MDRC}",
    xticks= -90:10:90,
    xgridvisible=false,
    ygridvisible=false
)

for fname in readdir(hyperbolic_dir)
    label = splitext(fname)[1]
    label == "1.0" && continue
    data = load_solver_data(joinpath(hyperbolic_dir, fname))
    mdrc_data = calc_mdrc(data)
    display(data.params)
    
    lines!(ax, mdrc_data.θs, mdrc_data.inc_p[:, 1], label=L"a = %$(label)")
end
axislegend(ax)
save("plots/hyberbolic_type1.pdf", fig)

hyperbolic2_dir = "output/sweep/hyperbolic_type2"

fig = Figure(fontsize=24, size=(800, 600))
ax = Axis(fig[1, 1],
    title = "MDRC",
    xlabel=L"\theta_s [deg]",
    ylabel=L"\text{MDRC}",
    xticks= -90:10:90,
    xgridvisible=false,
    ygridvisible=false
)

for fname in readdir(hyperbolic2_dir)
    label = splitext(fname)[1]
    label == "1.0" && continue
    data = load_solver_data(joinpath(hyperbolic2_dir, fname))
    mdrc_data = calc_mdrc(data)
    display(data.params)
    
    lines!(ax, mdrc_data.θs, mdrc_data.inc_p[:, 1], label=L"a = %$(label)")
end
axislegend(ax)
save("plots/hyberbolic_type2.pdf", fig)