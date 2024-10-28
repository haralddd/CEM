
using CairoMakie
using LaTeXStrings

using CairoMakie
using LaTeXStrings
CairoMakie.activate!(type="svg")
function plot_mdrc(; coh, incoh, θ0s=[0, 10, 20], colors=[:red, :blue, :green], name="default")

    mkpath("plots")
    θs = range(-90, 90, length=size(incoh, 1))

    # Incoherent MDRC
    fig = Figure(fontsize=24)
    ax = Axis(fig[1, 1], xlabel=L"$\theta_s$ [deg]", ylabel=L"Incoh. MDRC $$")
    for i in axes(incoh, 2)
        θ0 = θ0s[i]
        y = incoh[:, i]
        lines!(ax,
            θs, y,
            label=L"$\theta_0=$%$θ0",
            linestyle=:solid,
            linewidth=1,
            color=colors[i])
        vlines!(ax,
            [-θ0, θ0],
            color=colors[i],
            linestyle=:dot,
            linewidth=1)
    end
    axislegend()
    save("plots/$(name)_incoh.pdf", fig)

    # Coherent MDRC
    fig = Figure(fontsize=24)
    ax = Axis(fig[1, 1], xlabel=L"$\theta_s$ [deg]", ylabel=L"Coh. MDRC $$")
    for i in axes(coh, 2)
        θ0 = θ0s[i]
        y = coh[:, i]
        lines!(ax,
            θs, y,
            label=L"$\theta_0=$%$θ0",
            linestyle=:solid,   
            linewidth=1,
            color=colors[i])    
        vlines!(ax,
            [-θ0, θ0],
            color=colors[i],
            linestyle=:dot,
            linewidth=1)
    end
    axislegend()
    save("plots/$(name)_coh.pdf", fig)

    println("Saved plots to \'plots/$name\'")
end
