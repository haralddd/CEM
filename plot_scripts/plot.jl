using RayleighSolver

using CairoMakie
using LaTeXStrings
using Colors
CairoMakie.activate!(type="svg")

const tol = 1.0e-10

function /(x::AbstractString, y::AbstractString)
    return joinpath(x, y)
end

const DEFAULT_PATH = joinpath(@__DIR__, "..")
const DEFAULT_OUTPUT = DEFAULT_PATH / "output"
const DEFAULT_PLOTDIR = DEFAULT_PATH / "plots"
mkpath(DEFAULT_PLOTDIR)
mkpath(DEFAULT_OUTPUT)

colors = [:darkgreen, :royalblue, :firebrick, :darkorange]
markers = [:circle, :rect, :diamond, :utriangle]
markersize = 16
general_fontsize = 24
label_fontsize = 32
axis_fontsize = 40
ylabelpad = 10.0
linewidth = 1.0

# Direct plotting is used instead of these helper functions
# for better control over polarization-specific visualization

function save_mdrc_plots(data, θ0s, θs, ylabel, file_prefix, folder)
    # Create stacked plots with shared x-axis
    n_plots = size(data, 2)
    fig = Figure(fontsize=general_fontsize, size=(600, 300 * n_plots))
    
    # Use a color scheme that's distinct and print-friendly
    global colors, markers
    if n_plots > 4
        colors = distinguishable_colors(n_plots, [colorant"white", colorant"black"], dropseed=true)
    end
    
    # Create axes with linkxaxes
    fig_axes = []
    
    for i in 1:n_plots
        is_bottom = i == n_plots  # Bottom plot in the stack
        
        # Create plot
        ax = Axis(fig[i, 1], 
            xlabel = is_bottom ? L"\theta_s\ \text{[deg]}" : "",
            ylabel = ylabel,
            ylabelpadding = ylabelpad,
            ylabelsize = axis_fontsize,
            xlabelsize = axis_fontsize,
            xticks = -90:30:90,
            xminorticks = -90:10:90,
            xminorticksvisible = true,
            xminorgridvisible = true,
            xminorgridstyle = (:dot, :dense),
            xminorgridcolor = (:black, 0.12),
            xticklabelsvisible = is_bottom,
        )
        
        if i > 1
            linkxaxes!(fig_axes[1], ax)
        end
        push!(fig_axes, ax)
        
        # Plot the data
        θ0 = θ0s[i]
        θ0_label = round(θ0, digits=1)
        
        # Plot vertical lines at the incident angle
        vlines!(ax, [-θ0, θ0], color = (:black, 0.5), linestyle = :dash, linewidth = 1.0)
        
        # Plot the data
        lines!(ax, θs, data[:, i], linewidth = linewidth, color = :red)

        text!(ax, 0.95, 0.5, text = L"\theta_0=%$θ0_label^{\circ}", 
              fontsize = general_fontsize, align = (:right, :center), space = :relative)
    end
    
    # Save the stacked plot
    save(folder / "$(file_prefix)_stacked.pdf", fig)
    
    # Also create a combined plot where all angles are on the same axis
    fig_combined = Figure(fontsize=general_fontsize, size=(800, 600))
    ax_combined = Axis(fig_combined[1, 1], 
        xlabel = L"\theta_s\ \text{[deg]}", 
        ylabel = ylabel, 
        ylabelsize = axis_fontsize,
        xlabelsize = axis_fontsize,
        xticks = -90:30:90,
        xminorticks = -90:10:90,
        xminorticksvisible = true,
        xminorgridvisible = true,
        xminorgridstyle = (:dot, :dense),
        xminorgridcolor = (:black, 0.12),
        ylabelpadding = ylabelpad,
    )
    
    comb_max = 0.0
    # Plot all angles in the same figure with a legend
    for i in axes(data, 2)
        ys = data[:, i]
        θ0 = θ0s[i]
        θ0_label = round(θ0, digits=1)
        global colors, markers, markersize
        ln = lines!(ax_combined, θs, ys, linewidth = linewidth, color = colors[i])
        
        # Plot vertical lines at incident angles
        # vlines!(ax_combined, [-θ0, θ0], color = (colors[i], 0.3), 
        #         linestyle = :dash, linewidth = 1.0)
        
        ymax, maxidx = findmax(ys)
        sc = scatterlines!(ax_combined, θs[maxidx], ymax, color=colors[i], marker=markers[i], markersize=markersize,
            label=L"\theta_0=%$θ0_label^{\circ}")
        translate!(sc, 0, 0, 1000-ymax)
        translate!(ln, 0, 0, 1000-ymax)
    end
    
    # Add legend
    axislegend(ax_combined, labelsize=label_fontsize)
    
    # Set y limits based on all data
    _min = minimum(data)
    _max = maximum(data)
    _max += 0.5 * abs(_max) + tol
    ylims!(ax_combined, (_min, _max))
    
    # Save the combined plot
    save(folder / "$(file_prefix)_combined.pdf", fig_combined)
end

function save_mdtc_plots_p(data, θtes, θs, θ0s, ylabel, file_prefix, folder)
    # P-polarization function for extraordinary waves
    
    # Create stacked plots with shared x-axis
    n_plots = size(data, 2)
    fig = Figure(fontsize=general_fontsize, size=(600, 300 * n_plots))
    
    # Use a color scheme that's distinct and print-friendly
    global colors, markers, markersize
    if n_plots > 4
        colors = distinguishable_colors(n_plots, [colorant"white", colorant"black"], dropseed=true)
    end
    
    # Create axes with linkxaxes
    fig_axes = []
    
    for i in 1:n_plots
        is_bottom = i == n_plots  # Bottom plot in the stack
        
        # Create plot
        ax = Axis(fig[i, 1], 
            xlabel = is_bottom ? L"\theta_t\ \text{[deg]}" : "",
            ylabel = ylabel,
            ylabelsize = axis_fontsize,
            xlabelsize = axis_fontsize,
            xticks = -90:30:90,
            xminorticks = -90:10:90,
            xminorticksvisible = true,
            xminorgridvisible = true,
            xminorgridstyle = (:dot, :dense),
            xminorgridcolor = (:black, 0.12),
            xticklabelsvisible = is_bottom,
            ylabelpadding = ylabelpad,
        )
        
        if i > 1
            linkxaxes!(fig_axes[1], ax)
        end
        push!(fig_axes, ax)
        
        # Get the incident and extraordinary transmission angles
        θte = θtes[i]
        θ0_label = round(i ≤ length(θ0s) ? θ0s[i] : θte, digits=1)  # Use incident angle if available
        θt_label = round(θte, digits=1)  # p-pol uses extraordinary
        

        # For p-polarization, show extraordinary direction
        vlines!(ax, [-θte, θte], color = (:black, 0.5), linestyle = :dash, linewidth = 1.0)
        
        # Plot the data
        lines!(ax, θs, data[:, i], linewidth = linewidth, color = :red)

        # Add text annotation with normalized coordinates (0,0) to (1,1)
        # text!(ax, 0.05, 0.9, text = L"\theta_0=%$θ0_label^{\circ}", 
        #       fontsize = general_fontsize, align = (:left, :top), space = :relative)
        text!(ax, 0.95, 0.5, text = L"\theta_\mathrm{te}=%$θt_label^{\circ}", 
              fontsize = general_fontsize, align = (:right, :center), space = :relative)
        
    end
    
    # Save the stacked plot
    save(folder / "$(file_prefix)_stacked.pdf", fig)
    
    # Also create a combined plot where all angles are on the same axis
    fig_combined = Figure(fontsize=general_fontsize, size=(800, 600))
    ax_combined = Axis(fig_combined[1, 1], 
        xlabel = L"\theta_t\ \text{[deg]}", 
        ylabel = ylabel, 
        ylabelsize = axis_fontsize,
        xlabelsize = axis_fontsize,
        xticks = -90:30:90,
        xminorticks = -90:10:90,
        xminorticksvisible = true,
        xminorgridvisible = true,
        xminorgridstyle = (:dot, :dense),
        xminorgridcolor = (:black, 0.12),
        ylabelpadding = ylabelpad,
    )
    
    # Plot all angles in the same figure with a legend
    for i in axes(data, 2)
        ys = data[:, i]
        θte = θtes[i]
        θ0_label = round(i ≤ length(θ0s) ? θ0s[i] : θte, digits=1)
        
        # Use extraordinary angle for p-polarization
        θt_label = round(θte, digits=1)
        
        global colors, markers, markersize
        ln = lines!(ax_combined, θs, ys, linewidth = linewidth, color = colors[i])
        
        # For p-polarization, show extraordinary direction
        # vlines!(ax_combined, [-θte, θte], color = (colors[i], 0.3), 
        #         linestyle = :dash, linewidth = 1.0)
        maxidx = argmax(abs.(ys))
        ymax = abs(ys[maxidx])
        sc = scatterlines!(ax_combined, θs[maxidx], ys[maxidx], color=colors[i], marker=markers[i], markersize=markersize,
            label=L"\theta_0=%$θ0_label^{\circ},\ \theta_\mathrm{te}=%$θt_label^{\circ}")
        translate!(sc, 0, 0, 1000-ymax)
        translate!(ln, 0, 0, 1000-ymax)
    end
    
    # Add legend
    axislegend(ax_combined, labelsize=label_fontsize)
    
    # Set y limits based on all data
    _min = minimum(data)
    _max = maximum(data)
    _max += 0.5 * abs(_max) + tol
    ylims!(ax_combined, (_min, _max))
    
    # Save the combined plot
    save(folder / "$(file_prefix)_combined.pdf", fig_combined)
end

function save_mdtc_plots_s(data, θtos, θs, θ0s, ylabel, file_prefix, folder)
    # S-polarization function for ordinary waves
    
    # Create stacked plots with shared x-axis
    n_plots = size(data, 2)
    fig = Figure(fontsize=general_fontsize, size=(600, 300 * n_plots))
    
    # Use a color scheme that's distinct and print-friendly
    global colors, markers, markersize
    if n_plots > 4
        colors = distinguishable_colors(n_plots, [colorant"white", colorant"black"], dropseed=true)
    end
    
    # Create axes with linkxaxes
    fig_axes = []
    
    for i in 1:n_plots
        is_bottom = i == n_plots  # Bottom plot in the stack
        
        # Create plot
        ax = Axis(fig[i, 1], 
            xlabel = is_bottom ? L"\theta_t\ \text{[deg]}" : "",
            ylabel = ylabel,
            ylabelsize = axis_fontsize,
            xlabelsize = axis_fontsize,
            xticks = -90:30:90,
            xminorticks = -90:10:90,
            xminorticksvisible = true,
            xminorgridvisible = true,
            xminorgridstyle = (:dot, :dense),
            xminorgridcolor = (:black, 0.12),
            xticklabelsvisible = is_bottom,
            ylabelpadding = ylabelpad,
        )
        
        if i > 1
            linkxaxes!(fig_axes[1], ax)
        end
        push!(fig_axes, ax)
        
        # Get the incident and ordinary transmission angles
        θto = θtos[i]
        θ0_label = round(i ≤ length(θ0s) ? θ0s[i] : θto, digits=1)  # Use incident angle if available
        θt_label = round(θto, digits=1)  # s-pol uses ordinary
        

        # For s-polarization, show ordinary direction
        vlines!(ax, [-θto, θto], color = (:black, 0.5), linestyle = :dash, linewidth = 1.0)
        
        # Plot the data
        lines!(ax, θs, data[:, i], linewidth = linewidth, color = :red)

        # Add text annotation with normalized coordinates (0,0) to (1,1)
        # text!(ax, 0.05, 0.9, text = L"\theta_0=%$θ0_label^{\circ}", 
        #       fontsize = general_fontsize, align = (:left, :top), space = :relative)
        text!(ax, 0.95, 0.5, text = L"\theta_\mathrm{to}=%$θt_label^{\circ}", 
              fontsize = general_fontsize, align = (:right, :center), space = :relative)
        
    end
    
    # Save the stacked plot
    save(folder / "$(file_prefix)_stacked.pdf", fig)
    
    # Also create a combined plot where all angles are on the same axis
    fig_combined = Figure(fontsize=general_fontsize, size=(800, 600))
    ax_combined = Axis(fig_combined[1, 1], 
        xlabel = L"\theta_t\ \text{[deg]}", 
        ylabel = ylabel, 
        ylabelsize = axis_fontsize,
        xlabelsize = axis_fontsize,
        xticks = -90:30:90,
        xminorticks = -90:10:90,
        xminorticksvisible = true,
        xminorgridvisible = true,
        xminorgridstyle = (:dot, :dense),
        xminorgridcolor = (:black, 0.12),
        ylabelpadding = ylabelpad,
    )
    
    # Plot all angles in the same figure with a legend
    for i in axes(data, 2)
        ys = data[:, i]
        θto = θtos[i]
        θ0_label = round(i ≤ length(θ0s) ? θ0s[i] : θto, digits=1)
        
        # Use ordinary angle for s-polarization
        θt_label = round(θto, digits=1)
        ln = lines!(ax_combined, θs, ys, linewidth = linewidth, color = colors[i])
        
        # For s-polarization, show ordinary direction
        # vlines!(ax_combined, [-θto, θto], color = (colors[i], 0.3), 
        #         linestyle = :dash, linewidth = 1.0)
        
        maxidx = argmax(abs.(ys))
        ymax = abs(ys[maxidx])
        sc = scatterlines!(ax_combined, θs[maxidx], ys[maxidx], color=colors[i], marker=markers[i], markersize=markersize,
            label=L"\theta_0=%$θ0_label^{\circ},\ \theta_\mathrm{to}=%$θt_label^{\circ}")
        translate!(sc, 0, 0, 1000-ymax)
        translate!(ln, 0, 0, 1000-ymax)
    end
    
    # Add legend
    axislegend(ax_combined, labelsize=label_fontsize)
    
    # Set y limits based on all data
    _min = minimum(data)
    _max = maximum(data)
    _max += 0.5 * abs(_max) + tol
    ylims!(ax_combined, (_min, _max))
    
    # Save the combined plot
    save(folder / "$(file_prefix)_combined.pdf", fig_combined)
end

function make_plots(data::SolverData, fname="default", dir="plots")

    folder = dir / fname
    mkpath(folder)

    gen_ylabel(rt_str, pol_str, ang_str, type_str) = L"\langle\partial %$(rt_str)_{%$pol_str}/\partial\theta_{%$ang_str}\rangle_{\mathrm{%$type_str}}"
    
    mdrc = calc_mdrc(data)
    save_mdrc_plots(mdrc.coh_p, mdrc.θ0s, mdrc.θs, gen_ylabel("R", "p", "s", "coh"), "mdrc-p-coh", folder)
    save_mdrc_plots(mdrc.inc_p, mdrc.θ0s, mdrc.θs, gen_ylabel("R", "p", "s", "incoh"), "mdrc-p-incoh", folder)
    save_mdrc_plots(mdrc.coh_s, mdrc.θ0s, mdrc.θs, gen_ylabel("R", "s", "s", "coh"), "mdrc-s-coh", folder)
    save_mdrc_plots(mdrc.inc_s, mdrc.θ0s, mdrc.θs, gen_ylabel("R", "s", "s", "incoh"), "mdrc-s-incoh", folder)

    if data.solver_type == :full
        mdtc = calc_mdtc(data)
        save_mdtc_plots_p(mdtc.coh_p, mdtc.θtes, mdtc.θtps, mdrc.θ0s, gen_ylabel("T", "p", "t", "coh"), "mdtc-p-coh", folder)
        save_mdtc_plots_p(mdtc.inc_p, mdtc.θtes, mdtc.θtps, mdrc.θ0s, gen_ylabel("T", "p", "t", "incoh"), "mdtc-p-incoh", folder)
        save_mdtc_plots_s(mdtc.coh_s, mdtc.θtos, mdtc.θtss, mdrc.θ0s, gen_ylabel("T", "s", "t", "coh"), "mdtc-s-coh", folder)
        save_mdtc_plots_s(mdtc.inc_s, mdtc.θtos, mdtc.θtss, mdrc.θ0s, gen_ylabel("T", "s", "t", "incoh"), "mdtc-s-incoh", folder)
    end

    @info "Saved plots to $(folder)"
end

function cli_plot_main(arg)
    filename = ""
    try filename = files[parse(Int64, arg)]
    catch
        filename = arg
    end

    if isfile(DEFAULT_OUTPUT / filename)
        data = load_solver_data(DEFAULT_OUTPUT / filename)
    else
        @info "File '$filename' not found in $(DEFAULT_OUTPUT). Looking for matches."
        files = readdir(DEFAULT_OUTPUT)
        for file in files
            if occursin(arg, file)
                data = load_solver_data(DEFAULT_OUTPUT / file)
                @info "Using file '$file'"
                break
            end
        end
    end

    display(data.params)

    make_plots(data, filename, DEFAULT_PLOTDIR)
end

if (abspath(PROGRAM_FILE) == @__FILE__)
    files = readdir(DEFAULT_OUTPUT)
    if length(ARGS) < 1
        println("Usage: julia plot.jl <filename || index>")
        println("List of output files in $(DEFAULT_OUTPUT):")

        for (idx, file) in enumerate(files)
            println("$idx: $file")
        end
    else
        cli_plot_main(ARGS[1])
    end
end