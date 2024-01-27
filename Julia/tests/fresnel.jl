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

    Δθ = 0.5
    θs = 0.0:Δθ:90.0-Δθ
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
        # Draw new surface sample
        generate!(sp_p.ys, rp_p, surf_t)
        generate!(sp_s.ys, rp_s, surf_t)

        # Calculate the invariant part of Npk (depends on k)
        N_invariant!(Npk_p_pre, rp_p, ks[i])
        N_invariant!(Npk_s_pre, rp_s, ks[i])

        # Solve for R
        solve!(sp_p, rp_p, Mpk_p_pre, Npk_p_pre, kis[i])
        solve!(sp_s, rp_s, Mpk_s_pre, Npk_s_pre, kis[i])

        res_p[i] = (sp_p.Npk .|> abs |> maximum) ./ 2π
        res_s[i] = (sp_s.Npk .|> abs |> maximum) ./ 2π
    end

    fig1 = Figure(; size=(500, 400))
    ax1 = Axis(fig1[1, 1], xlabel=L"Reflected angle $θ$", ylabel=L"Reflection coefficient $R$")
    scatter!(ax1, θs, res_p .^ 2, label="p", color=:blue, marker=:circle)
    scatter!(ax1, θs, res_s .^ 2, label="s", color=:red, marker=:rect)
    lines!(ax1, θs, Rp.(θs, ε), label="p (Fresnel)", linestyle=:solid, color=(:blue, 0.5))
    lines!(ax1, θs, Rs.(θs, ε), label="s (Fresnel)", linestyle=:solid, color=(:red, 0.5))
    axislegend(ax1; position=:lt)
    tightlimits!(ax1)

    # Plot error between 1 and sum of Fresnel coefficients
    fig2 = Figure(; size=(500, 400))
    ax2 = Axis(fig2[1, 1], xlabel=L"Reflected angle $θ$", ylabel=L"Error, $|1 - R_p + R_s|$")

    lines!(ax2, θs, abs.(1.0 .- (res_p .^ 2 .+ res_s .^ 2)))
    tightlimits!(ax2)

    return fig1, fig2
end

fig1, fig2 = test_fresnel(ε=2.25)
display(fig1)
display(fig2)
save("plots/fresnel_glass.pdf", fig1)
save("plots/fresnel_glass_error.pdf", fig2)
fig1, fig2 = test_fresnel(ε=-17.5 + 0.48im)
display(fig1)
display(fig2)
save("plots/fresnel_silver.pdf", fig1)
save("plots/fresnel_silver_error.pdf", fig2)