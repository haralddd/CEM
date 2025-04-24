using RayleighSolver
using CairoMakie
using LaTeXStrings

const eps_perp = -7.5 + 0.24im
const surf_rms = 10.0e-9
const surf_corlen = 100.0e-9
const θs = [10.0]
const λ = 457.9e-9
const Nx = 4096
const Lx = 200
const Ni = 10
const ens_iters = 10000

surface = GaussianSurface(surf_rms, surf_corlen)
above = Vacuum()

# Parameter ranges to sweep

metal_ratios = [0.5, 1.0, 1.5]
mkpath("output/sweep/metallic")
mdrcs = []
for ar in metal_ratios
    below = Uniaxial(eps_perp, eps_perp*ar, 1.0 + 0.0im, 1.0 + 0.0im)
    params = Parameters(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)

    data = SolverData(params, ens_iters, :reduced)
    solve_ensemble!(data)
    mdrc = calc_mdrc(data)
    push!(mdrcs, mdrc)
end
fig = Figure()
ax = fig[1, 1] = Axis(fig)

for (ar, mdrc) in zip(metal_ratios, mdrcs)
    lines!(ax, mdrc.θs, mdrc.inc_p[:], label=L"\mathrm{Inc.\ MDRC} (\nu = p,\ a = %$ar)")
end
axislegend(ax)
display(fig)
save("output/sweep/metallic.pdf", fig)



hyperbolic_ratios = [-0.5, -1.0, -1.5]
mkpath("output/sweep/hyperbolic")
mdrcs = []
for ar in hyperbolic_ratios
    below = Uniaxial(eps_perp, eps_perp*ar, 1.0 + 0.0im, 1.0 + 0.0im)
    params = Parameters(surf=surface, above=above, below=below, Nx=Nx, Lx=Lx, lambda=λ, θs=θs)

    data = SolverData(params, ens_iters, :reduced)
    solve_ensemble!(data)
    mdrc = calc_mdrc(data)
    push!(mdrcs, mdrc)
end
fig = Figure()
ax = fig[1, 1] = Axis(fig)

for (ar, mdrc) in zip(hyperbolic_ratios, mdrcs)
    lines!(ax, mdrc.θs, mdrc.inc_p[:], label=L"\mathrm{Inc.\ MDRC} (u = p, a = %$ar)")
end
axislegend(ax)
display(fig)
save("output/sweep/hyperbolic.pdf", fig)