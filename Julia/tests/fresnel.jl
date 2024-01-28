push!(LOAD_PATH, "Julia/RayleighSolver/")
using RayleighSolver
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


function test_fresnel(; input="silver")

    rp_p, sp_p, generator_p! = load_solver_config("input/$(input)_p_fresnel.txt")
    rp_s, sp_s, generator_s! = load_solver_config("input/$(input)_s_fresnel.txt")

    # Calc the invariant part of Mpk
    display("Calculating invariant parts of Mpk")
    Mpk_p_pre = Array{ComplexF64,3}(undef, length(rp_p.ps), length(rp_p.qs), rp_p.Ni + 1)
    Npk_p_pre = Array{ComplexF64,3}(undef, length(rp_p.ps), length(rp_p.kis), rp_p.Ni + 1)

    Mpk_s_pre = Array{ComplexF64,3}(undef, length(rp_s.ps), length(rp_s.qs), rp_s.Ni + 1)
    Npk_s_pre = Array{ComplexF64,3}(undef, length(rp_s.ps), length(rp_s.kis), rp_s.Ni + 1)

    @time M_invariant!(Mpk_p_pre, rp_p)
    @time M_invariant!(Mpk_s_pre, rp_s)

    # Calculate the invariant part of Npk (depends on k)
    @time N_invariant!(Npk_p_pre, rp_p)
    @time N_invariant!(Npk_s_pre, rp_s)

    # Check for undefined behaviour
    @assert all(isfinite.(Mpk_p_pre))
    @assert all(isfinite.(Mpk_s_pre))
    @assert all(isfinite.(Npk_p_pre))
    @assert all(isfinite.(Npk_s_pre))

    # Generate surface (flat)
    generator_p!(sp_p.ys)
    generator_s!(sp_s.ys)

    # Solve for R
    @time solve!(sp_p, rp_p, Mpk_p_pre, Npk_p_pre)
    @time solve!(sp_s, rp_s, Mpk_s_pre, Npk_s_pre)

    res_p = [sp_p.Npk[:, i] .|> abs2 |> maximum for i in eachindex(rp_p.kis)]
    res_s = [sp_s.Npk[:, i] .|> abs2 |> maximum for i in eachindex(rp_s.kis)]

    return res_p, res_s, asind.(rp_p.qs[rp_p.kis])
end


function plot_fresnel(; res_p, res_s, θs, εμ, title="silver")

    fig1 = Figure(; size=(500, 400))
    ax1 = Axis(fig1[1, 1], xlabel=L"Reflected angle $θ$", ylabel=L"Reflection coefficient $R$")
    scatter!(ax1, θs, res_p, label="p (Numerical)", color=:blue, marker=:circle)
    scatter!(ax1, θs, res_s, label="s (Numerical)", color=:red, marker=:rect)
    lines!(ax1, θs, Rp.(θs, εμ), label="p (Analytical)", linestyle=:solid, color=:blue)
    lines!(ax1, θs, Rs.(θs, εμ), label="s (Analytical)", linestyle=:solid, color=:red)
    axislegend(ax1; position=:lt)
    limits!(ax1, 0, 90, 0.97, 1.001)

    # tightlimits!(ax1)

    display(fig1)

    # Plot error between 1 and sum of Fresnel coefficients
    fig2 = Figure(; size=(500, 400))
    ax2 = Axis(fig2[1, 1], xlabel=L"Reflected angle $θ$", ylabel=L"Absolute error$$")

    lines!(ax2, θs, abs.(res_p .- Rp.(θs, εμ)), label="p", linestyle=:solid, color=:blue)
    lines!(ax2, θs, abs.(res_s .- Rs.(θs, εμ)), label="s", linestyle=:solid, color=:red)
    tightlimits!(ax2)
    display(fig2)
    save("plots/fresnel_$(title).pdf", fig1)
    save("plots/fresnel_$(title)_error.pdf", fig2)
end

res_p, res_s, θs = test_fresnel()
ε = -17.5 + 0.48im

plot_fresnel(; res_p=res_p, res_s=res_s, θs=θs, εμ=ε * 1.0, title="silver")