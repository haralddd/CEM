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
function custom_plot!(ax, xs, ys, x0, color, shape)
    lines!(ax,
        xs, ys,
        label=L"$\theta_0=$%$x0",
        linestyle=:solid,
        linewidth=1,
        color=color)
    vlines!(ax,
        [-x0, x0],
        color=color,
        linestyle=:dot,
        linewidth=1)
    maxidx = argmax(ys)
    scatter!(ax,
        [xs[maxidx]],
        [ys[maxidx]],
        color=color,
        marker=shape,
        markersize=16,
    )
    return
end

function plot_mdrc(data::SolverData, fname="default", dir="plots")

    colors = [:red, :green, :blue, :purple, :black, :orange, :cyan, :magenta]
    shapes = [:circle, :cross, :rect, :diamond, :utriangle, :dtriangle, :ltriangle, :rtriangle]

    coh = data.out.coherent
    incoh = data.out.incoherent
    θ0s = asind.(data.spa.ks)
    θs = asind.(data.out.qs)


    # Incoherent MDRC
    fig = Figure(fontsize=24)
    ax = Axis(fig[1, 1], xlabel=L"$\theta_s$ [deg]", ylabel=L"Incoh. MDRC $$")
    for i in axes(incoh, 2)
        θ0 = θ0s[i]
        ys = incoh[:, i]
        custom_plot!(ax, θs, ys, θ0, colors[i], shapes[i])
    end
    axislegend()
    save(joinpath(dir, fname * "_incoh.pdf"), fig)

    # Coherent MDRC
    fig = Figure(fontsize=24)
    ax = Axis(fig[1, 1], xlabel=L"$\theta_s$ [deg]", ylabel=L"Coh. MDRC $$")
    for i in axes(coh, 2)
        θ0 = θ0s[i]
        ys = coh[:, i]
        custom_plot!(ax, θs, ys, θ0, colors[i], shapes[i])
    end
    axislegend()
    save(joinpath(dir, fname * "_coh.pdf"), fig)

    @info "Saved plots to $(DEFAULT_PLOTDIR / fname)"
end

function cli_plot(arg)
    filename = ""
    arg = ARGS[1]
    try filename = files[parse(Int64, arg)]
    catch
        filename = arg
    end

    plot_mdrc(load_solver_data(DEFAULT_OUTPUT / filename), filename, DEFAULT_PLOTDIR)
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
        cli_plot(ARGS[1])
    end
end