module RayleighSolver

using LinearAlgebra
using FFTW
using Base.Threads
LinearAlgebra.BLAS.set_num_threads(1)

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
    solve!, run_threaded



end # module RayleighSolver
