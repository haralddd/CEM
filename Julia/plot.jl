push!(LOAD_PATH, "$(@__DIR__)/RayleighSolver/")
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
function custom_plot!(ax, xs, ys, x0)
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

function save_plots(data, θ0s, θs, label, folder)
    for i in axes(data, 2)
        fig = Figure(fontsize=24)
        ax = Axis(fig[1, 1], xlabel=L"$\theta_s$ [deg]", ylabel=L"\langle \partial_{\theta_s} R\rangle_{\mathrm{%$label}}")
        θ0 = θ0s[i]
        ys = data[:, i]
        custom_plot!(ax, θs, ys, θ0)
        axislegend()
        save(folder / "$(label)_$(i).pdf", fig)
    end
end

function calc_mdrc_and_save_plots(data::SolverData, fname="default", dir="plots")

    folder = dir / fname
    mkpath(folder)

    ks = data.spa.ks

    qs, mdrc = get_qs_and_mdrc(data)
    # full_qs = data.spa.qs
    # σ² = data.out.σ²
    # κ = data.out.κ

    theta0s = asind.(ks)
    thetas = asind.(qs)

    # P-polarization
    # Incoherent MDRC
    save_plots(mdrc.coh_p, theta0s, thetas, "p,coh", folder)
    save_plots(mdrc.inc_p, theta0s, thetas, "p,incoh", folder)
    save_plots(mdrc.coh_s, theta0s, thetas, "s,coh", folder)
    save_plots(mdrc.inc_s, theta0s, thetas, "s,incoh", folder)

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
        @info "File '$filename' not found in $(DEFAULT_OUTPUT). Look for matches."
        files = readdir(DEFAULT_OUTPUT)
        for file in files
            if occursin(arg, file)
                data = load_solver_data(DEFAULT_OUTPUT / file)
                @info "Using file '$file'"
                break
            end
        end
    end

    calc_mdrc_and_save_plots(data, filename, DEFAULT_PLOTDIR)
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