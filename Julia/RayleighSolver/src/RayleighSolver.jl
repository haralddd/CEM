module RayleighSolver

using LinearAlgebra
using FFTW

include("setup.jl")
include("solver.jl")

export
    # Setup functions
    Polarization, p, s,
    polarization_from_string,
    RayleighParams,
    params_as_string,
    SurfPreAlloc,
    SurfType, flat, gaussian, singlebump,
    generate!,


    # Solver functions
    α, α0,
    M_ker, N_ker,
    M_invariant!, N_invariant!,
    solve!



end # module RayleighSolver
