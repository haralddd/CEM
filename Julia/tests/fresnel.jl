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


function test_fresnel(; ε=2.25)
    surf_t::SurfType = flat::SurfType

    Δθ = 1.0
    θs = 0.0:Δθ:90.0
    Nq = 2^10

    rp_p = RayleighParams(
        ν=p,
        Nq=Nq,
        ε=ε
    )

    rp_s = RayleighParams(
        ν=s,
        Nq=Nq,
        ε=ε
    )

    # Calc the invariant part of Mpk
    display("Calculating invariant parts of Mpk")
    Mpk_p_pre = Array{ComplexF64,3}(undef, length(rp_p.ps), length(rp_p.qs), rp_p.Ni + 1)
    Npk_p_pre = Matrix{ComplexF64}(undef, length(rp_p.ps), rp_p.Ni + 1)

    Mpk_s_pre = Array{ComplexF64,3}(undef, length(rp_s.ps), length(rp_s.qs), rp_s.Ni + 1)
    Npk_s_pre = Matrix{ComplexF64}(undef, length(rp_s.ps), rp_s.Ni + 1)

    @time M_invariant!(Mpk_p_pre, rp_p)
    @time M_invariant!(Mpk_s_pre, rp_s)

    # Check for undefined behaviour
    @assert all(isfinite.(Mpk_p_pre))
    @assert all(isfinite.(Mpk_s_pre))

    kis = [searchsortedfirst(rp_p.qs, sind(θs[i]), rev=true) for i in eachindex(θs)] |> unique
    ks = rp_p.qs[kis]
    θs = asind.(ks)

    sp_p = SurfPreAlloc(rp_p, surf_t)
    sp_s = SurfPreAlloc(rp_s, surf_t)

    # Results in Fresnel coefficients
    res_p = Vector{Float64}(undef, length(ks))
    res_s = Vector{Float64}(undef, length(ks))

    @time for i in eachindex(ks)

        # Calculate the invariant part of Npk (depends on k)
        N_invariant!(Npk_p_pre, rp_p, ks[i])
        N_invariant!(Npk_s_pre, rp_s, ks[i])

        # Solve for R
        solve!(sp_p, rp_p, Mpk_p_pre, Npk_p_pre, kis[i])
        solve!(sp_s, rp_s, Mpk_s_pre, Npk_s_pre, kis[i])

        res_p[i] = (sp_p.Npk .|> abs2 |> maximum)
        res_s[i] = (sp_s.Npk .|> abs2 |> maximum)
    end

    return res_p, res_s, θs
end


function plot_fresnel(; ε=2.25, title="glass")

    res_p, res_s, θs = test_fresnel(ε=ε)

    fig1 = Figure(; size=(500, 400))
    ax1 = Axis(fig1[1, 1], xlabel=L"Reflected angle $θ$", ylabel=L"Reflection coefficient $R$")
    scatter!(ax1, θs, res_p, label="p (Numerical)", color=:blue, marker=:circle)
    scatter!(ax1, θs, res_s, label="s (Numerical)", color=:red, marker=:rect)
    lines!(ax1, θs, Rp.(θs, ε), label="p (Analytical)", linestyle=:solid, color=:blue)
    lines!(ax1, θs, Rs.(θs, ε), label="s (Analytical)", linestyle=:solid, color=:red)
    axislegend(ax1; position=:lt)
    limits!(ax1, 0, 90, 0, 1.05)

    # tightlimits!(ax1)

    display(fig1)

    # Plot error between 1 and sum of Fresnel coefficients
    fig2 = Figure(; size=(500, 400))
    ax2 = Axis(fig2[1, 1], xlabel=L"Reflected angle $θ$", ylabel=L"Absolute error$$")

    lines!(ax2, θs, abs.(res_p .- Rp.(θs, ε)), label="p", linestyle=:solid, color=:blue)
    lines!(ax2, θs, abs.(res_s .- Rs.(θs, ε)), label="s", linestyle=:solid, color=:red)
    tightlimits!(ax2)
    display(fig2)
    save("plots/fresnel_$(title).pdf", fig1)
    save("plots/fresnel_$(title)_error.pdf", fig2)
end

plot_fresnel(ε=2.25, title="glass")
plot_fresnel(ε=-17.5 + 0.48im, title="silver")