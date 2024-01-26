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

    @time pre_M_invariant!(Mpk_p_pre, rp_p)
    @time pre_M_invariant!(Mpk_s_pre, rp_s)

    # Check for undefined behaviour
    @assert all(isfinite.(Mpk_p_pre))
    @assert all(isfinite.(Mpk_s_pre))

    kis = [searchsortedfirst(rp_p.qs, sind(θs[i]), rev=true) for i in eachindex(θs)] |> unique
    ks = rp_p.qs[kis]
    θs = asind.(ks)

    # Results in Fresnel coefficients
    rs_p = Vector{Float64}(undef, length(ks))
    rs_s = Vector{Float64}(undef, length(ks))

    @time for i in eachindex(ks)
        sp_p = SurfPreAlloc(rp_p, surf_t)
        sp_s = SurfPreAlloc(rp_s, surf_t)

        pre_N_invariant!(Npk_p_pre, rp_p, ks[i])
        pre_N_invariant!(Npk_s_pre, rp_s, ks[i])

        # Check for undefined behaviour
        Np_check = isfinite.(Npk_p_pre)
        Ns_check = isfinite.(Npk_s_pre)
        if !all(Np_check) || !all(Ns_check)
            idxs_p = findall(x -> !x, Np_check)
            idxs_s = findall(x -> !x, Ns_check)
            display("Npk_p_pre ($idxs_p) or Npk_s_pre ($idxs_s) is not finite")

            display("$(Npk_p_pre[idxs_p]), $(Npk_s_pre[idxs_s])")
            display("Iteration: $i")

            error("Undefined behaviour")
        end

        # Solve and insert the specular reflection coefficient
        solve_pre!(sp_p, rp_p, Mpk_p_pre, Npk_p_pre, kis[i])
        solve_pre!(sp_s, rp_s, Mpk_s_pre, Npk_s_pre, kis[i])

        rs_p[i] = sp_p.R .|> abs2 |> maximum
        rs_s[i] = sp_s.R .|> abs2 |> maximum
    end


    display("Plotting Fresnel coefficients")
    fig1 = Figure(; size=(500, 400))
    ax1 = Axis(fig1[1, 1], xlabel=L"Reflected angle $θ$", ylabel=L"Reflection coefficient $R$")
    scatter!(ax1, θs, rs_p, label="p", color=:blue, marker=:circle)
    scatter!(ax1, θs, rs_s, label="s", color=:red, marker=:rect)
    lines!(ax1, θs, Rp.(θs, ε), label="p (Fresnel)", linestyle=:solid, color=(:blue, 0.5))
    lines!(ax1, θs, Rs.(θs, ε), label="s (Fresnel)", linestyle=:solid, color=(:red, 0.5))
    axislegend(ax1; position=:lt)
    tightlimits!(ax1)

    # Plot error between 1 and sum of Fresnel coefficients
    fig2 = Figure(; size=(500, 400))
    ax2 = Axis(fig2[1, 1], xlabel=L"Reflected angle $θ$", ylabel=L"Error, $|1 - R_p + R_s|$")

    lines!(ax2, θs, abs.(1.0 .- (rs_p .^ 2 .+ rs_s .^ 2)))
    tightlimits!(ax2)

    return fig1, fig2
end

fig1, fig2 = test_fresnel(ε=2.25 + 1e-6im)
display(fig1)
display(fig2)
save("plots/fresnel_glass.pdf", fig1)
save("plots/fresnel_glass_error.pdf", fig2)
fig1, fig2 = test_fresnel(ε=-17.5 + 0.48im)
display(fig1)
display(fig2)
save("plots/fresnel_silver.pdf", fig1)
save("plots/fresnel_silver_error.pdf", fig2)

