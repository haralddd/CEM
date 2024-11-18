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
import Base.parse
import Base.show
import Base.display
import Base.convert

include("Material.jl")
export Material, Vacuum, Isotropic, UniaxialCrystal
export _A, ptilde, alpha, alpha0

include("RandomSurface.jl")
export RandomSurface, FlatSurface, GaussianSurface
export SingleBumpSurface, RectangularSurface
export scale

include("SimParams.jl")
export Polarization, PolarizationP, PolarizationS
export SimParams
export get_angles, get_scale, get_scaled_params

include("SimPrealloc.jl")
export SystemPreAlloc, SimPrealloc
export precompute!, validate

include("random_surface_generator.jl")
export generate_surface!

include("solver.jl")
export SimOutput, SolverData, DataMDRC
export solve_single!, solve_MDRC!, precompute!, observe, observe!
export calc_mdrc

include("utils.jl")
export show, display, parse, convert
export save_spa_config, load_spa_config
export save_solver_data, load_solver_data
export save_ensemble_iters, load_ensemble_iters
export surface_prompt, config_creation_prompt, default_config_creation
export default_params_for_surface_testing

end # module RayleighSolver
