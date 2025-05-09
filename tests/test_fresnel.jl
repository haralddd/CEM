include("testconfig.jl")
using Statistics
using CairoMakie
using LaTeXStrings

function Rs(θ, εμ)
    # Fresnel reflection coefficient
    n1 = 1 # relative vacuum
    n2 = sqrt(εμ)
    abs2((n1 * cosd(θ) - n2 * sqrt(1 - (n1 / n2 * sind(θ))^2)) /
         (n1 * cosd(θ) + n2 * sqrt(1 - (n1 / n2 * sind(θ))^2)))
end

function Rp(θ, εμ)
    # Fresnel reflection coefficient
    n1 = 1 # relative vacuum
    n2 = sqrt(εμ)
    abs2((n1 * sqrt(1 - (n1 / n2 * sind(θ))^2) - n2 * cosd(θ)) /
         (n1 * sqrt(1 - (n1 / n2 * sind(θ))^2) + n2 * cosd(θ)))
end

function plot_fresnel(; res_p, res_s, θs, εμ, title="silver")

    fig1 = Figure(; size=(500, 400))
    ax1 = Axis(fig1[1, 1], xlabel=L"Reflected angle $θ$", ylabel=L"Reflection coefficient $R$")
    scatter!(ax1, θs, res_p, label="p (Numerical)", color=:blue, marker=:circle)
    scatter!(ax1, θs, res_s, label="s (Numerical)", color=:red, marker=:rect)
    lines!(ax1, θs, Rp.(θs, εμ), label="p (Analytical)", linestyle=:solid, color=:blue)
    lines!(ax1, θs, Rs.(θs, εμ), label="s (Analytical)", linestyle=:solid, color=:red)
    axislegend(ax1; position=:lt)
    # limits!(ax1, 0, 90, 0.96, 1.001)

    display(fig1)

    # Plot error between 1 and sum of Fresnel coefficients
    fig2 = Figure(; size=(500, 400))
    ax2 = Axis(fig2[1, 1], xlabel=L"Reflected angle $θ$", ylabel=L"Absolute error$$")

    lines!(ax2, θs, abs.(res_p .- Rp.(θs, εμ)), label="p", linestyle=:solid, color=:blue)
    lines!(ax2, θs, abs.(res_s .- Rs.(θs, εμ)), label="s", linestyle=:solid, color=:red)
    tightlimits!(ax2)
    mkpath("plots")
    save("plots/fresnel_$(title).pdf", fig1)
    save("plots/fresnel_$(title)_error.pdf", fig2)
end

function test_fresnel(;type=:silver)
    if type == :silver
        data = config_fresnel_silver(1000)
    elseif type == :glass
        data = config_fresnel_glass(1000)
    end

    prealloc = Preallocated(data.params)
    precompute!(data.precomputed, data.params)
    generate_surface!(prealloc, data.params)
    solve_single!(prealloc, data)

    pNpk = prealloc.PNpk
    sNpk = prealloc.SNpk

    res_p = [pNpk[:, i] .|> abs2 |> maximum for i in axes(pNpk, 2)]
    res_s = [sNpk[:, i] .|> abs2 |> maximum for i in axes(sNpk, 2)]

    εμ = data.params.below.eps * data.params.below.mu

    plot_fresnel(; res_p=res_p, res_s=res_s, θs=data.params.θs, εμ=εμ, title=string(type))

    return nothing
end

test_fresnel(type=:silver)
test_fresnel(type=:glass)