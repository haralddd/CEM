using RayleighSolver
using CairoMakie
using LaTeXStrings
using DataFrames
# Use the analytical expression for a flat surface to plot the reflection coefficients
# Compare with the numerical results with a flat surface

function R(a, a0, κpa)
    pre = 1 / (a/κpa + a0)
    return pre * (-a/κpa + a0)
end

function T(a, a0, κpa)
    pre = 1 / (a/κpa + a0)
    return pre * (2a0)
end

function ema_eps_para(f, eps_m, eps_d)
    return f * eps_m + (1 - f)*eps_d
end

function ema_eps_perp(f, eps_m, eps_d)
    return eps_m*eps_d / ((1-f)*eps_m + f*eps_d)
end


@info "Flat surface test"
θs = 0.0:1.0:90.0
surf = FlatSurface()
above = Vacuum()

# TiO2
eps_tio2_para = 3.62+0.0im
eps_tio2_perp = 6.84+0.0im

# Gold and glass layered
# # eps_gold = -11.740+1.2611im
# # eps_glass = 2.25+0.0im
# # eps_gg_perp = ema_eps_perp(0.165, eps_gold, eps_glass)
# # eps_gg_para = ema_eps_para(0.165, eps_gold, eps_glass)


# Silver and glass layered
eps_silver = -17.5+0.48im
eps_glass = 2.25 + 0.0im
eps_sg_perp = ema_eps_perp(0.2, eps_silver, eps_glass)
eps_sg_para = ema_eps_para(0.2, eps_silver, eps_glass)
below = Uniaxial(eps_sg_perp, eps_sg_para, 1.0+0.0im, 1.0+0.0im)


Nx = 2*2048
data = SolverData(Parameters(surf=surf, θs=θs, above=above, below=below, Nx=Nx, Ni=1))
θs2 = unique(data.params.θs)
data = SolverData(Parameters(surf=surf, θs=θs2, above=above, below=below, Nx=Nx, Ni=1))
prealloc = Preallocated(data.params)
precompute!(data.precomputed, data.params)
generate_surface!(prealloc, data.params)
solve_single!(prealloc, data)


ks = sind.(θs)
ap = alpha_p.(ks, A(below), eps_para)
as = alpha_s.(ks, eps_para)
a0 = alpha0.(ks)

Rp = R.(ap, a0, below.eps_para)
Rs = R.(as, a0, below.mu_para)

plt = plot(θs, abs2.(Rp), label=L"Analysis EMA, $\nu=p$", 
    legendfontsize=10, color=:blue, ylabel=L"R^\nu", xlabel=L"\theta_i")
plot!(θs, abs2.(Rs), label=L"Analysis EMA, $\nu=s$", color=:red)

kis = data.params.kis
scatter!(θs2, [abs2.(prealloc.PNpk[ki, i] for (i,ki) in enumerate(kis))], label=L"RRE Solver EMA, $\nu=p$", color=:blue, markersize=2)
scatter!(θs2, [abs2.(prealloc.SNpk[ki, i] for (i, ki) in enumerate(kis))], label=L"RRE Solver EMA, $\nu=s$", color=:red, markersize=2)

# PyLlama values
pyllama_Rs = CSV.read("output/PyLlama-RSs_a100nm.csv", DataFrame, header=0, delim=',', types=Float64)
pyllama_Rp = CSV.read("output/PyLlama-RPs_a100nm.csv", DataFrame, header=0, delim=',', types=Float64)

pl_θs = rad2deg.(pyllama_Rp[:, 1])
pl_Rp = pyllama_Rp[:, 2]
pl_Rs = pyllama_Rs[:, 2]
plot!(pl_θs, pl_Rp, label=L"PyLlama $\nu=p$", color=:blue, linestyle=:dot)
plot!(pl_θs, pl_Rs, label=L"PyLlama $\nu=s$", color=:red, linestyle=:dot)
# ylims!(0.0, 1.05)
xticks!(0:10:90)

savefig(plt, "plots/fresnel_ema_a100nm.pdf")

pl_ks = sind.(pl_θs)
ap = alpha_p.(pl_ks, A(below), eps_para)
as = alpha_s.(pl_ks, eps_para)
a0 = alpha0.(pl_ks)
Rp = R.(ap, a0, below.eps_para)
Rs = R.(as, a0, below.mu_para)

ΔRp = abs.(abs2.(Rp) .- pl_Rp)
ΔRs = abs.(abs2.(Rs) .- pl_Rs)

plt_err = plot(pl_θs, ΔRp, label=L"$\nu = p$", color=:blue, legend=:topleft, ylabel=L"|\Delta R^\nu|", xlabel=L"\theta_i")
plot!(pl_θs, ΔRs, label=L"$\nu = s$", color=:red)
# plot!([1], [NaN], label=L"$\left|\frac{\Delta R^p}{R^p}\right|$", color=:blue, linestyle=:dot)
# plot!([1], [NaN], label=L"$\left|\frac{\Delta R^s}{R^s}\right|$", color=:red, linestyle=:dot)

# ΔRp_rel = 100*ΔRp ./ pl_Rp
# ΔRs_rel = 100*ΔRs ./ pl_Rs

# ax = twinx()
# ylabel!(ax, "%")
# plot!(ax, pl_θs, ΔRp_rel, label=nothing, color=nothing, linestyle=:dot)
# plot!(ax, pl_θs, ΔRs_rel, label=nothing, color=nothing, linestyle=:dot)

savefig(plt_err, "plots/fresnel_ema_error_a100nm.pdf")

