using RayleighSolver
using CairoMakie
using LaTeXStrings
using DataFrames, CSV
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
# eps_tio2_para = 3.62+0.0im
# eps_tio2_perp = 6.84 + 0.0im
# below = Uniaxial(eps_tio2_perp, eps_tio2_para, 1.0 + 0.0im, 1.0 + 0.0im)

# Gold and glass layered
# # eps_gold = -11.740+1.2611im
# # eps_glass = 2.25+0.0im
# # eps_gg_perp = ema_eps_perp(0.165, eps_gold, eps_glass)
# # eps_gg_para = ema_eps_para(0.165, eps_gold, eps_glass)


# Silver and glass layered
eps_silver = -17.5+0.48im
eps_glass = 2.25 + 0.0im
eps_sg_perp = ema_eps_perp(0.5, eps_silver, eps_glass)
eps_sg_para = ema_eps_para(0.5, eps_silver, eps_glass)
below = Uniaxial(eps_sg_perp, eps_sg_para, 1.0+0.0im, 1.0+0.0im)


Nx = 2*2048
paramconf = ParametersConfig(surf=surf, θs=θs, above=above, below=below, Nx=Nx, Ni=1)
params = Parameters(paramconf)
θs = unique(params.θs)
paramconf = ParametersConfig(surf=surf, θs=θs, above=above, below=below, Nx=Nx, Ni=1)

data = SolverData(paramconf, 1, :reduced)
params = data.params
prealloc = Preallocated(data)
precomputed = Precomputed(data)
precompute!(precomputed, data)
generate_surface!(prealloc, data.params)
solve_single!(prealloc, precomputed, data)
observe!(data, prealloc, 1)


ks = sind.(θs)
ap = alpha_p.(ks, Ref(below))
as = alpha_s.(ks, Ref(below))
a0 = alpha0.(ks)

Rp = R.(ap, a0, below.eps_para)
Rs = R.(as, a0, below.mu_para)


fig = Figure(size=(800, 600), fontsize=36)
ax = fig[1, 1] = Axis(fig, 
    xlabel=L"\theta_i\ [^\circ]", ylabel=L"R_f^\nu",
    xticks=(0:30:90))
lines!(ax, θs, abs2.(Rp), label=L"Analysis, $\nu=p$", color=:blue)
lines!(ax, θs, abs2.(Rs), label=L"Analysis, $\nu=s$", color=:red)

kis = data.params.kis
rp = [data.Rp.A²[ki,i] for (i,ki) in enumerate(kis)]
rs = [data.Rs.A²[ki,i] for (i,ki) in enumerate(kis)]
scatter!(ax, θs, rp, label=L"Solver, $\nu=p$", color=:blue)
scatter!(ax, θs, rs, label=L"Solver, $\nu=s$", color=:red)
axislegend(ax,position=:lt)
display(fig)
save("plots/fresnel_scatter_sg.pdf", fig)

# PyLlama values
pyllama_Rs_10 = CSV.read("output/PyLlama-RSs_a10nm.csv", DataFrame, header=0, delim=',', types=Float64)
pyllama_Rp_10 = CSV.read("output/PyLlama-RPs_a10nm.csv", DataFrame, header=0, delim=',', types=Float64)
pyllama_Rs_100 = CSV.read("output/PyLlama-RSs_a100nm.csv", DataFrame, header=0, delim=',', types=Float64)
pyllama_Rp_100 = CSV.read("output/PyLlama-RPs_a100nm.csv", DataFrame, header=0, delim=',', types=Float64)

pl_θs = rad2deg.(pyllama_Rp_10[:, 1])
pl_Rp_10 = pyllama_Rp_10[:, 2]
pl_Rs_10 = pyllama_Rs_10[:, 2]

pl_θs = rad2deg.(pyllama_Rp_100[:, 1])
pl_Rp_100 = pyllama_Rp_100[:, 2]
pl_Rs_100 = pyllama_Rs_100[:, 2]

fig = Figure(size=(800, 600), fontsize=24)
ax = fig[1, 1] = Axis(fig, 
    xlabel=L"\theta_i\ [^\circ]", ylabel=L"R_f^\nu",
    xticks=(0:30:90))
lines!(ax, θs, abs2.(Rp), label=L"EMA, $\nu=p$", color=:blue)
lines!(ax, θs, abs2.(Rs), label=L"EMA, $\nu=s$", color=:red)
lines!(ax, pl_θs, pl_Rp_10, label=L"PyLlama, $a=10$nm, $\nu=p$", color=:blue, linestyle=:dot)
lines!(ax, pl_θs, pl_Rs_10, label=L"PyLlama, $a=10$nm, $\nu=s$", color=:red, linestyle=:dot)
lines!(ax, pl_θs, pl_Rp_100, label=L"PyLlama, $a=100$nm, $\nu=p$", color=:blue, linestyle=:dash)
lines!(ax, pl_θs, pl_Rs_100, label=L"PyLlama, $a=100$nm, $\nu=s$", color=:red, linestyle=:dash)
ylims!(ax, 0.96, nothing)

axislegend(ax, position=:lb)
display(fig)
save("plots/fresnel_ema_vs_pl.pdf", fig)

