"""Implements optical scattering under the reduced Rayleigh equations
Effectively solves Maxwell's equations for a singularly polarized wave
on a partially symmetric rough surface. Satisfying the boundary conditions

The reduced Rayleigh equations are a set of coupled integral equations
which assume that far field scattering conditions (singular direction, up/down in 1D)
can be used all the way down to the rough surface boundary even though for strongly
rough surfaces one can get multiple scattering events.
"""
module RayleighSolver

using LinearAlgebra
using FFTW
using Statistics
using Random: Xoshiro, randn!
using JLD2, FileIO
using ProgressBars
using Roots
import Base.parse
import Base.show
import Base.display
import Base.convert
import Base.+

include("Material.jl")
export Material, Vacuum, Isotropic, Uniaxial
export A, alpha, alpha0, alpha_p, alpha_s

include("RandomSurface.jl")
export RandomSurface, FlatSurface, GaussianSurface
export SingleBumpSurface, RectangularSurface, scale

include("Parameters.jl")
export Polarization, PolarizationP, PolarizationS
export Parameters
export get_angles, get_scale, get_scaled_params

include("Precomputed.jl")
export Precomputed, precompute!, validate

include("Preallocated.jl")
export SolverData, Preallocated, Results, get_A, get_A², observe, observe!

include("random_surface_generator.jl")
export generate_surface!


include("compute_t.jl")
export compute_t_matrix, compute_t_matrix_full

include("solver_reduced_isotropic.jl")
include("solver_reduced_uniaxial.jl")
include("solver_full_uniaxial.jl")
export PlotData
export solve_single_full!, solve_single_reduced!, solve_single!, precompute!

include("energy.jl")
export energy_conservation, energy_ratio

include("solver_ensemble.jl")
export solve_ensemble!

include("calc_plotdata.jl")
export MdrcPlotData, MdtcPlotData, calc_mdrc, calc_mdtc


include("utils.jl")
export show, display, parse, convert
export save_spa_config, load_spa_config
export save_solver_data, load_solver_data
export save_ensemble_iters, load_ensemble_iters
export surface_prompt, config_creation_prompt, default_config_creation
export default_params_for_surface_testing

end # module RayleighSolver
