using RayleighSolver

using CairoMakie
using LaTeXStrings
CairoMakie.activate!(type="svg")

function /(x::AbstractString, y::AbstractString)
    return joinpath(x, y)
end

const DEFAULT_PATH = "$(@__DIR__)/"
const DEFAULT_OUTPUT = DEFAULT_PATH / "output"
const DEFAULT_PLOTDIR = DEFAULT_PATH / "plots"
mkpath(DEFAULT_PLOTDIR)
mkpath(DEFAULT_OUTPUT)
function mdrc_plot!(ax, xs, ys, x0)
    x0 = round(x0, digits=1)
    lines!(ax,
        xs, ys,
        linestyle=:solid,
        linewidth=0.5,
        color=:red)
    vlines!(ax,
        [-x0, x0],
        label=L"\theta_0=\pm %$x0",
        color=:red,
        linestyle=:dot,
        linewidth=1.0)
    return
end

function mdtc_plot!(ax, xs, ys, x0s)

    @assert length(x0s) == 2 "Should be exactly 2 Fresnel angles for mdtc plot"

    lines!(ax,
        xs, ys,
        linestyle=:solid,
        linewidth=0.5,
        color=:red)

    
    xo = round(x0s[1], digits=1) # Ordinary wave direction
    vlines!(ax,
        [-xo, xo],
        label=L"\theta_o=\pm %$xo",
        color=:red,
        linestyle=:dot,
        linewidth=1.0)

    xe = round(x0s[2], digits=1) # Extraordinary wave direction
    vlines!(ax,
        [-xe, xe],
        label=L"\theta_e=\pm %$xe",
        color=:red,
        linestyle=:dash,
        linewidth=1.0)

    return
end

function save_mdrc_plots(data, θ0s, θs, title, file_prefix, folder)
    for i in axes(data, 2)
        fig = Figure(fontsize=24)
        ax = Axis(fig[1, 1], title = title, xlabel=L"$\theta_s$ [deg]", ylabel=L"\text{MDRC}")
        θ0 = θ0s[i]
        ys = data[:, i]
        ylims!(ax, (0,maximum(ys) .* 1.2))
        mdrc_plot!(ax, θs, ys, θ0)
        axislegend()
        save(folder / "$(file_prefix)_$(i).pdf", fig)
    end
end

function save_mdtc_plots(data, θ0s, θs, title, file_prefix, folder)
    for i in axes(data, 2)
        fig = Figure(fontsize=24)
        ax = Axis(fig[1, 1], title = title, xlabel=L"$\theta_t$ [deg]", ylabel=L"\text{MDTC}")
        θ0 = θ0s[i]
        ys = data[:, i]
        ylims!(ax, (0,maximum(ys) .* 1.2))
        mdtc_plot!(ax, θs, ys, θ0)
        axislegend()
        save(folder / "$(file_prefix)_$(i).pdf", fig)
    end
end

function make_plots(data::SolverData, fname="default", dir="plots")

    folder = dir / fname
    mkpath(folder)

    mdrc_data = calc_mdrc(data)
    θ0s = mdrc_data.θ0s
    θss = mdrc_data.θs

    cp = mdrc_data.coh_p
    ip = mdrc_data.inc_p
    cs = mdrc_data.coh_s
    is = mdrc_data.inc_s

    @debug "θ0s: $θ0s"
    @debug "θss: $θss"

    # P-polarization
    # Incoherent MDRC
    save_mdrc_plots(cp, θ0s, θss, L"\text{Coherent MDRC, }\nu = p", "mdrc-p-coh", folder)
    save_mdrc_plots(ip, θ0s, θss, L"\text{Incoherent MDRC, }\nu = p", "mdrc-p-incoh", folder)
    save_mdrc_plots(cs, θ0s, θss, L"\text{Coherent MDRC, }\nu = s", "mdrc-s-coh", folder)
    save_mdrc_plots(is, θ0s, θss, L"\text{Incoherent MDRC, }\nu = s", "mdrc-s-incoh", folder)

    @info "Saved plots to $(folder)"
end

function make_plots(data::SolverData{Parameters{_S,Vacuum,Uniaxial}}, fname="default", dir="plots") where _S

    folder = dir / fname
    mkpath(folder)

    mdrc = calc_mdrc(data)
    mdtc = calc_mdtc(data)
    θss = mdrc.θ0s
    θts = mdtc.θ0s
    θs = mdrc.θs

    @debug "θs: $θs"
    @debug "θss: $θss"
    @debug "θts: $θts"

    dq = data.params.dq
    @info "∑MDRC_s = $((sum(mdrc.coh_s) + sum(mdrc.inc_s))*dq)"
    @info "∑MDRC_p = $((sum(mdrc.coh_p) + sum(mdrc.inc_p))*dq)" 
    @info "∑MDTC_s = $((sum(mdtc.coh_s) + sum(mdtc.inc_s))*dq)"
    @info "∑MDTC_p = $((sum(mdtc.coh_p) + sum(mdtc.inc_p))*dq)"

    # P-polarization
    # Incoherent MDRC
    save_mdrc_plots(mdrc.coh_p, θss, θs, L"\text{Coherent MDRC, }\nu = p", "mdrc-p-coh", folder)
    save_mdrc_plots(mdrc.inc_p, θss, θs, L"\text{Incoherent MDRC, }\nu = p", "mdrc-p-incoh", folder)
    save_mdrc_plots(mdrc.coh_s, θss, θs, L"\text{Coherent MDRC, }\nu = s", "mdrc-s-coh", folder)
    save_mdrc_plots(mdrc.inc_s, θss, θs, L"\text{Incoherent MDRC, }\nu = s", "mdrc-s-incoh", folder)

    save_mdtc_plots(mdtc.coh_p, θts, θs, L"\text{Coherent MDTC, }\nu = p", "mdtc-p-coh", folder)
    save_mdtc_plots(mdtc.inc_p, θts, θs, L"\text{Incoherent MDTC, }\nu = p", "mdtc-p-incoh", folder)
    save_mdtc_plots(mdtc.coh_s, θts, θs, L"\text{Coherent MDTC, }\nu = s", "mdtc-s-coh", folder)
    save_mdtc_plots(mdtc.inc_s, θts, θs, L"\text{Incoherent MDTC, }\nu = s", "mdtc-s-incoh", folder)


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