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
export A, ptilde, alpha, alpha0, alpha_p, alpha_s

include("RandomSurface.jl")
export RandomSurface, FlatSurface, GaussianSurface
export SingleBumpSurface, RectangularSurface
export scale

include("Parameters.jl")
export Polarization, PolarizationP, PolarizationS
export Parameters
export get_angles, get_scale, get_scaled_params

include("Precomputed.jl")
export Precomputed
export precompute!, validate

include("Preallocated.jl")
export Preallocated, Results, get_R, get_T, get_R², get_T²

include("random_surface_generator.jl")
export generate_surface!

include("SolverData.jl")
export SolverData

include("solver_reduced.jl")
include("solver_full.jl")
export PlotData
export solve_single!, precompute!, observe, observe!

include("energy.jl")
export energy_conservation, energy_ratio

include("solver_ensemble.jl")
export solve_MDRC!, calc_mdrc, calc_mdtc


include("utils.jl")
export show, display, parse, convert
export save_spa_config, load_spa_config
export save_solver_data, load_solver_data
export save_ensemble_iters, load_ensemble_iters
export surface_prompt, config_creation_prompt, default_config_creation
export default_params_for_surface_testing

end # module RayleighSolver
