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
function mdrc_plot!(ax, xs, ys, θ0)
    θ0 = round(θ0, digits=2)
    lines!(ax,
        xs, ys,
        linestyle=:solid,
        linewidth=0.5,
        color=:red)
    vlines!(ax,
        [-θ0, θ0],
        label=L"\theta_0=\pm %$θ0",
        color=:red,
        linestyle=:dot,
        linewidth=1.0)
    return
end

function mdtc_plot!(ax, xs, ys, θto, θte)

    lines!(ax,
        xs, ys,
        linestyle=:solid,
        linewidth=0.5,
        color=:red)

    
    θo = round(θto, digits=2) # Ordinary wave direction
    vlines!(ax,
        [-θo, θo],
        label=L"\theta_o=\pm %$θo",
        color=:red,
        linestyle=:dot,
        linewidth=1.0)

    θe = round(θte, digits=2) # Extraordinary wave direction
    vlines!(ax,
        [-θe, θe],
        label=L"\theta_e=\pm %$θe",
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
        ylims!(ax, (minimum(ys), maximum(ys) .* 1.3))
        mdrc_plot!(ax, θs, ys, θ0)
        axislegend()
        save(folder / "$(file_prefix)_$(i).pdf", fig)
    end
end

function save_mdtc_plots(data, θtos, θtes, θs, title, file_prefix, folder)
    for i in axes(data, 2)
        fig = Figure(fontsize=24)
        ax = Axis(fig[1, 1], title = title, xlabel=L"$\theta_t$ [deg]", ylabel=L"\text{MDTC}")
        θto = θtos[i]
        θte = θtes[i]
        ys = data[:, i]
        ylims!(ax, (minimum(ys), maximum(ys) .* 1.3))
        mdtc_plot!(ax, θs, ys, θto, θte)
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
    mdtc_p, mdtc_s = calc_mdtc(data)

    @debug "mdrc: $mdrc"
    @debug "mdtc_p: $mdtc_p"
    @debug "mdtc_s: $mdtc_s"

    @info "∑MDRC_s = $(sum(mdrc.coh_s) + sum(mdrc.inc_s))"
    @info "∑MDRC_p = $(sum(mdrc.coh_p) + sum(mdrc.inc_p))"
    @info "∑MDTC_s = $(sum(mdtc_s.coh) + sum(mdtc_s.inc))"
    @info "∑MDTC_p = $(sum(mdtc_p.coh) + sum(mdtc_p.inc))"

    @info "mdtc_p.θtos: $(mdtc_p.θtos)"
    @info "mdtc_p.θtes: $(mdtc_p.θtes)"

    @info "mdtc_p.θs: $(mdtc_p.θs[1]):$(mdtc_p.θs[end])"
    @info "mdtc_s.θs: $(mdtc_s.θs[1]):$(mdtc_s.θs[end])"

    @info "qs len $(length(data.params.qs))"
    @info "θtp len $(length(mdtc_p.θs))"
    @info "θts len $(length(mdtc_s.θs))"
    

    # P-polarization
    # Incoherent MDRC
    save_mdrc_plots(mdrc.coh_p, mdrc.θ0s, mdrc.θs, L"\text{Coherent MDRC, }\nu = p", "mdrc-p-coh", folder)
    save_mdrc_plots(mdrc.inc_p, mdrc.θ0s, mdrc.θs, L"\text{Incoherent MDRC, }\nu = p", "mdrc-p-incoh", folder)
    save_mdrc_plots(mdrc.coh_s, mdrc.θ0s, mdrc.θs, L"\text{Coherent MDRC, }\nu = s", "mdrc-s-coh", folder)
    save_mdrc_plots(mdrc.inc_s, mdrc.θ0s, mdrc.θs, L"\text{Incoherent MDRC, }\nu = s", "mdrc-s-incoh", folder)

    save_mdtc_plots(mdtc_p.coh, mdtc_p.θtos, mdtc_p.θtes, mdtc_p.θs, L"\text{Coherent MDTC, }\nu = p", "mdtc-p-coh", folder)
    save_mdtc_plots(mdtc_p.inc, mdtc_p.θtos, mdtc_p.θtes, mdtc_p.θs, L"\text{Incoherent MDTC, }\nu = p", "mdtc-p-incoh", folder)
    save_mdtc_plots(mdtc_s.coh, mdtc_s.θtos, mdtc_s.θtes, mdtc_s.θs, L"\text{Coherent MDTC, }\nu = s", "mdtc-s-coh", folder)
    save_mdtc_plots(mdtc_s.inc, mdtc_s.θtos, mdtc_s.θtes, mdtc_s.θs, L"\text{Incoherent MDTC, }\nu = s", "mdtc-s-incoh", folder)


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