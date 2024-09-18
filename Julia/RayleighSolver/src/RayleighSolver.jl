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
import Base.parse
import Base.show
import Base.display
import Base.convert

include("simulation_prealloc.jl")
include("surface.jl")
include("setup.jl")
include("surface_generator.jl")
include("utils.jl")
include("solver.jl")
# From simulation_prealloc.jl
export SimulationPreAlloc

# From surface.jl
export SurfaceParams, FlatSurfaceParams, GaussianSurfaceParams, SingleBumpSurfaceParams, RectSurfaceParams
export scale

# From setup.jl
export Polarization, p, s
export Scale, real, natural
export RayleighParams, c0, get_angles, scaled_params, get_scale

# From surface_generator.jl
export generate!, correlation

# From solver.jl
export solve!, solve_MDRC!

# From utils.jl
export config_creation_prompt, default_config_creation, default_params_for_surface_testing
export save_to, load_rp_struct, load_rp_desc
export parse, show, display, convert # Overloaded base methods for most objects

end # module RayleighSolver
