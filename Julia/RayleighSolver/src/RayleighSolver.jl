module RayleighSolver

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


    # Solver functions
    α, α0,
    M_ker, N_ker,
    M_invariant!, N_invariant!,
    M_invariant, N_invariant,
    pre_M_invariant!, pre_N_invariant!,
    solve_pre!, solve!



end # module RayleighSolver
