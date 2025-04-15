using RayleighSolver
using LinearAlgebra
using Plots
using LaTeXStrings
include("../src/tmatrix.jl")

"""
Test the T-matrix numerical integration implementation against the expected behavior
for a simple surface configuration.
"""

# Setup test parameters similar to test_RT_relation.jl
θ0 = [10.0]
material = Uniaxial(2.25 + 0.0im, 2.25 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im)
surf = GaussianSurface(10.0e-9, 100.0e-9)
above = Vacuum()

# Setup solver data
Nx = 2*1024  # Number of spatial points
data = SolverData(Parameters(surf=surf, θs=θ0, above=above, below=material, Nx=Nx, Ni=1))
prealloc = Preallocated(data.params)
precompute!(data.precomputed, data.params)
generate_surface!(prealloc, data.params)

# Solve for R-matrix
solve_single!(prealloc, data)
observe!(data.P_res, prealloc.PNpk, 1)
observe!(data.S_res, prealloc.SNpk, 1)

# Extract parameters needed for T-matrix calculation
qs = data.params.qs
ks = [data.params.k0 * sind(θ) for θ in θ0]
A = prealloc.A
μεpa = material.εxx / material.μxx

# Get R-matrix data
Rp = prealloc.PNpk
Rs = prealloc.SNpk

# Calculate T-matrix for p-polarization
println("Computing T-matrix for p-polarization...")
L = data.params.L
dp = 2π / L
p_values = [i * dp - π/dp for i in 0:(Nx-1)]

# For a single q,k pair as a test
q_idx = div(Nx, 2) + 1  # Middle q value
k_idx = 1  # First k value (from θ0)
q = qs[q_idx]
k = ks[k_idx]

# Extract R values for this q
R_values = Rp[q_idx, :]

# Compute single T-matrix element
T_element = compute_t_matrix(q, k, R_values, :p, A, μεpa, p_values, dp)
println("T^p($q|$k) = $T_element")

# Compute full T-matrix (can be time consuming for large Nx)
println("Computing full T-matrix (this may take a while for large Nx)...")
T_full_p = compute_t_matrix_full(qs, ks, Rp, :p, A, μεpa, Nx)
println("T-matrix computation complete.")

# Visualization
plot_q_idx = 1:10:Nx  # Subsample for plotting
heatmap(real.(T_full_p[plot_q_idx, 1]), 
        title="Real part of T-matrix for p-polarization", 
        xlabel="q index", ylabel="k index",
        color=:viridis)
savefig("t_matrix_real_p.png")

heatmap(imag.(T_full_p[plot_q_idx, 1]), 
        title="Imaginary part of T-matrix for p-polarization",
        xlabel="q index", ylabel="k index", 
        color=:viridis)
savefig("t_matrix_imag_p.png")

# Basic verification of energy conservation
# For energy conservation, we expect that |R|^2 + |T|^2 = 1 for lossless materials
R_energy = abs2(Rp[q_idx, k_idx])
T_energy = abs2(T_element)
energy_sum = R_energy + T_energy
println("Energy check for q=$q, k=$k:")
println("  |R|^2 = $R_energy")
println("  |T|^2 = $T_energy")
println("  |R|^2 + |T|^2 = $energy_sum (should be close to 1.0 for lossless systems)")
