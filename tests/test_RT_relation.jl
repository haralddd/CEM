using RayleighSolver
using Plots
using LaTeXStrings
using LinearAlgebra
using Statistics

# Test to verify the relation between transmission and reflection modes
# T^ν(q|k) = R^ν(q|k)e^{i(α(q,ω)+α_0(q,ω))ζ(x_1)} + 2πδ(q - k)e^{i(α_0(q,ω)-α_0(q,ω))ζ(x_1)}

@info "Testing relation between transmission and reflection modes"

δ(q) = q≈0.0 ? 1.0 : 0.0

# Function to calculate transmission modes from reflection modes
# based on the relation in the image
function calculate_T_from_R(R_qk, q, k, α_q, α0_q, ζ_x1)
    # First term: R^ν(q|k)e^{i(α(q,ω)+α_0(q,ω))ζ(x_1)}
    first_term = R_qk * exp(1im * (α_q + α0_q) * ζ_x1)
    
    # Second term: 2πδ(q - k)e^{i(α_0(q,ω)-α_0(q,ω))ζ(x_1)}
    # Note: The second exponential term simplifies to 1 since α_0(q,ω)-α_0(q,ω) = 0
    # The delta function is 1 when q = k, and 0 otherwise
    second_term = δ(q - k) * exp(1im * (α_q - α0_q) * ζ_x1)
    
    return first_term + second_term
end

# Setup test parameters
θ0 = [10.0]
material = Uniaxial(2.25 + 0.0im, 2.25 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im)
name = "glass"
surf = GaussianSurface(10.0e-9, 100.0e-9)
above = Vacuum()  # Medium above is vacuum

@info "Testing with $name below the interface"

# Setup solver data
Nx = 2*1024  # Number of spatial points
data = SolverData(Parameters(surf=surf, θs=θ0, above=above, below=material, Nx=Nx, Ni=1))
prealloc = Preallocated(data.params)
precompute!(data.precomputed, data.params)
generate_surface!(prealloc, data.params)
solve_single!(prealloc, data)
observe!(data.P_res, prealloc.PNpk, 1)
observe!(data.S_res, prealloc.SNpk, 1)

# Get the reflection coefficients from the solver
Rp_values = get_R(data.P_res)
Rs_values = get_R(data.S_res)

# Get the transmission coefficients from the solver
Tp_values = get_T(data.P_res)
Ts_values = get_T(data.S_res)

apqs = alpha_p.(data.params.qs, A(data.params.below), data.params.below.eps_para)
asqs = alpha_s.(data.params.qs, data.params.below.eps_para)

α0qs = alpha0.(data.params.qs)
ζ_x1 = 0.0

calculated_Tp_values = calculate_T_from_R.(Rp_values, data.params.qs, data.params.ks, apqs, α0qs, ζ_x1)
calculated_Ts_values = calculate_T_from_R.(Rs_values, data.params.qs, data.params.ks, asqs, α0qs, ζ_x1)

# Compare calculated T with actual T from solver
Tp_diff = calculated_Tp_values .- Tp_values
Ts_diff = calculated_Ts_values .- Ts_values

abs_err_p = abs.(Tp_diff)
abs_err_s = abs.(Ts_diff)

@info "Maximum absolute error for p-polarization: $(maximum(abs_err_p))"
@info "Maximum absolute error for s-polarization: $(maximum(abs_err_s))"

# Plot the results
plt_p = plot(abs.(Tp_values), label=L" \mathrm{Solver}\ T^p", 
    legendfontsize=10, color=:blue, ylabel=L"|T^p|", xlabel=L"q")
plot!(plt_p, abs.(calculated_Tp_values), label=L"\mathrm{Calculated}\ T^p", 
    linestyle=:dash, color=:red)
savefig(plt_p, "plots/RT_relation_$(name)_p.pdf")

plt_s = plot(abs.(Ts_values), label=L"\mathrm{Solver}\ T^s", 
    legendfontsize=10, color=:blue, ylabel=L"|T^s|", xlabel=L"q")
plot!(plt_s, abs.(calculated_Ts_values), label=L"\mathrm{Calculated}\ T^s", 
    linestyle=:dash, color=:red)
savefig(plt_s, "plots/RT_relation_$(name)_s.pdf")

# Plot the absolute difference
plt_diff = plot(abs_err_p, label=L"|T^p - T^p_{\text{calc}}|", 
    legendfontsize=10, color=:blue, ylabel="Absolute Error", xlabel=L"q", yscale=:log10)
plot!(plt_diff, abs_err_s, label=L"|T^s - T^s_{\text{calc}}|", 
    color=:red)
savefig(plt_diff, "plots/RT_relation_$(name)_abs_error.pdf")

# Calculate and plot the relative error
Tp_rel_error = abs_err_p ./ abs.(Tp_values) .* 100  # in percentage
Ts_rel_error = abs_err_s ./ abs.(Ts_values) .* 100  # in percentage

@info "Maximum relative error for p-polarization: $(maximum(Tp_rel_error))%"
@info "Maximum relative error for s-polarization: $(maximum(Ts_rel_error))%"

plt_rel_error = plot(Tp_rel_error, label=L"\nu=p", 
    legendfontsize=10, color=:blue, ylabel="Relative Error (%)", xlabel=L"q")
plot!(plt_rel_error, Ts_rel_error, label=L"\nu=s", 
    color=:red)
savefig(plt_rel_error, "plots/RT_relation_$(name)_rel_error.pdf")

