#=
    Includes all relevant definitions and includes for
    running a simulation
=#

const ε_AIR = 1.0005898


# Surface ensemble generation utils
include("SurfaceGen.jl")
using .SurfaceGen
export SurfaceEnsemble, Surface, mean_slope

# Source matrix generation
include("SourceMatrices.jl")
using .SourceMatrices
export create_A!, create_B!, SourceParams