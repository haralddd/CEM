using RayleighSolver
using CairoMakie
using LaTeXStrings
using LinearAlgebra
using Statistics
using FFTW

# Test to verify the relation between transmission and reflection modes
# T^ν(q|k) = R^ν(q|k)e^{i(α(q,ω)+α_0(q,ω))ζ(x_1)} + 2πδ(q - k)e^{i(α_0(q,ω)-α_0(q,ω))ζ(x_1)}

@info "Testing relation between transmission and reflection modes"

δ(q) = q ≈ 0.0 ? 1.0 : 0.0


function calc_Tp_Ts(Rp, Rs, pre::Preallocated, params::Parameters)
    qs = params.qs
    ps = params.ps
    ks = params.ks
    kis = params.kis

    ys = pre.ys
    Fys = pre.Fys
    sFys = pre.sFys
    FFT = params.FFT
    sFys_pqidxs = params.sFys_pqidxs

    Nq = size(qs, 1)
    Ni = params.Ni

    Ap = zeros(ComplexF64, Nq, Nq)
    As = zeros(ComplexF64, Nq, Nq)
    bp = zeros(ComplexF64, size(Rp))
    bs = zeros(ComplexF64, size(Rs))

    for n in 0:Ni
        for i in eachindex(Fys)
            Fys[i] = ys[i]^n
        end
        FFT * Fys
        fftshift!(sFys, Fys)

        for (kidx, kqidx) in enumerate(kis)
            a0 = alpha0(ks[kidx])
            for pidx in eachindex(ps)
                fourier_idx = sFys_pqidxs[pidx, kqidx]

                _bn = (-1.0im * a0)^n / factorial(n) * sFys[fourier_idx]
                bp[pidx, kidx] += _bn
                bs[pidx, kidx] += _bn
            end

            for qidx in eachindex(qs)
                ap = alpha_p(qs[qidx], params.below)
                as = alpha_s(qs[qidx], params.below)
                a0 = alpha0(qs[qidx])
                for pidx in eachindex(ps)
                    fourier_idx = sFys_pqidxs[pidx, qidx]

                    Ap[pidx, qidx] += (-1.0im * ap)^n * sFys[fourier_idx] / factorial(n)
                    As[pidx, qidx] += (-1.0im * as)^n * sFys[fourier_idx] / factorial(n)

                    bp[pidx, kidx] += (1.0im * a0)^n * sFys[fourier_idx] / factorial(n) * Rp[qidx, kidx]
                    bs[pidx, kidx] += (1.0im * a0)^n * sFys[fourier_idx] / factorial(n) * Rs[qidx, kidx]
                end
            end
        end
    end

    Tp = Ap \ bp
    Ts = As \ bs
    return Tp, Ts
end


# Setup test parameters
θ0 = [10.0]
material = Uniaxial(2.25 + 0.0im, 2.25 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im)
name = "glass"
surf = GaussianSurface(10.0e-9, 100.0e-9)
above = Vacuum()  # Medium above is vacuum

@info "Testing with $name below the interface"

# Setup solver data
Nx = 2 * 1024  # Number of spatial points
data = SolverData(Parameters(surf=surf, θs=θ0, above=above, below=material, Nx=Nx, Ni=1))
prealloc = Preallocated(data.params)
precomputed = Precomputed(data.params)
precompute!(precomputed, data.params)
generate_surface!(prealloc, data.params)
solve_single!(prealloc, precomputed, data)
half = size(prealloc.PNpk, 1) ÷ 2
Rp = prealloc.PNpk[1:half, :]
Rs = prealloc.SNpk[1:half, :]
Tp = prealloc.PNpk[half+1:end, :]
Ts = prealloc.SNpk[half+1:end, :]

Tp_calced, Ts_calced = calc_Tp_Ts(Rp, Rs, prealloc, data.params)

# Compare calculated T with actual T from solver
Tp_diff = Tp_calced .- Tp
Ts_diff = Ts_calced .- Ts

abs_err_p = abs.(Tp_diff)
abs_err_s = abs.(Ts_diff)

@info "Maximum absolute error for p-polarization: $(maximum(abs_err_p))"
@info "Maximum absolute error for s-polarization: $(maximum(abs_err_s))"

fig = Figure(size=(800, 600), fontsize=24)
ax = fig[1, 1] = Axis(fig, xlabel=L"q", ylabel=L"|T|", yscale=log10)
# Plot the results
lines!(ax, vec(abs.(Tp)), label=L" \mathrm{Solver}\ T^p", color=:blue)
lines!(ax, vec(abs.(Tp_calced)), label=L"\mathrm{Calculated}\ T^p",
    linestyle=:dash, color=:red)
axislegend(ax)
display(fig)
save("plots/RT_relation_$(name)_p.pdf", fig)

fig = Figure(size=(800, 600), fontsize=24)
ax = fig[1, 1] = Axis(fig, xlabel=L"q", ylabel=L"|T|", yscale=log10)
# Plot the results
lines!(ax, vec(abs.(Ts)), label=L" \mathrm{Solver}\ T^s", color=:blue)
lines!(ax, vec(abs.(Ts_calced)), label=L"\mathrm{Calculated}\ T^s",
    linestyle=:dash, color=:red)
axislegend(ax)
display(fig)
save("plots/RT_relation_$(name)_s.pdf", fig)

# Plot the absolute difference
fig = Figure(size=(800, 600), fontsize=24)
ax = fig[1, 1] = Axis(fig, xlabel=L"q", ylabel=L"|T - T_{\text{calc}}|", yscale=log10)
lines!(ax, vec(abs_err_p), label=L"|T^p - T^p_{\text{calc}}|", color=:blue)
lines!(ax, vec(abs_err_s), label=L"|T^s - T^s_{\text{calc}}|", color=:red)
axislegend(ax)
display(fig)
save("plots/RT_relation_$(name)_abs_error.pdf", fig)

# Calculate and plot the relative error
Tp_rel_error = abs_err_p ./ abs.(Tp) .* 100  # in percentage
Ts_rel_error = abs_err_s ./ abs.(Ts) .* 100  # in percentage

@info "Maximum relative error for p-polarization: $(maximum(Tp_rel_error))%"
@info "Maximum relative error for s-polarization: $(maximum(Ts_rel_error))%"

fig = Figure(size=(800, 600), fontsize=24)
ax = fig[1, 1] = Axis(fig, xlabel=L"q", ylabel=L"\mathrm{Relative\ error}\ [%]", yscale=log10)
lines!(ax, vec(Tp_rel_error), label=L"\nu=p", color=:blue)
lines!(ax, vec(Ts_rel_error), label=L"\nu=s", color=:red)
axislegend(ax)
display(fig)
save("plots/RT_relation_$(name)_rel_error.pdf", fig)

