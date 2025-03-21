using RayleighSolver
using Plots
using LaTeXStrings
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
eps_para = 3.62+0.0im
eps_perp = 6.84+0.0im

# Silver and glass layered
# eps_silver = -17.5+0.48im
# eps_glass = 2.25+0.0im
# eps_para = ema_eps_para(0.5, eps_silver, eps_glass)
# eps_perp = ema_eps_perp(0.5, eps_silver, eps_glass)
below = Uniaxial(eps_perp, eps_para, 1.0+0.0im, 1.0+0.0im)

data = SolverData(Parameters(surf=surf, θs=θs, above=above, below=below))
prealloc = Preallocated(data.params)
precompute!(data.precomputed, data.params)
generate_surface!(prealloc, data.params)
solve_single!(prealloc, data)


ks = sind.(θs)
ap = alpha_p.(ks, A(below), eps_perp)
as = alpha_s.(ks, eps_para)
a0 = alpha0.(ks)

Rp = R.(ap, a0, below.eps_para)
Rs = R.(as, a0, below.mu_para)

plt = plot(θs, abs2.(Rp), label=L"\textrm{Analytical}\ R^p", legendfontsize=12, color=:blue)
plot!(θs, abs2.(Rs), label=L"\textrm{Analytical}\ R^s", color=:red)
scatter!(θs, abs2.(Rp), label=L"\textrm{Numerical}\ R^p", color=:blue, markersize=2)
scatter!(θs, abs2.(Rs), label=L"\textrm{Numerical}\ R^s", color=:red, markersize=2)
# ylims!(0.0, 1.05)
xticks!(0:10:90)

savefig(plt, "plots/fresnel_tio2.pdf")
# savefig(plt, "plots/fresnel_layered_silver_glass.pdf")

display(plt)


    



