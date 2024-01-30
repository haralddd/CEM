module RayleighSolver

using LinearAlgebra
using FFTW
using Statistics
import Base.parse
import Base.string

include("setup.jl")
include("solver.jl")

export
    ## Setup
    # Enums and helper functions
    c0,
    Polarization, p, s,
    SurfType, flat, gaussian, singlebump, rect,
    polarization_from_string, surftype_from_string,
    Wg, gg, Wr, gr,

    # Structs used in computations
    Surface, RayleighParams,
    as_string, parse, SurfPreAlloc,
    flat_gen!, single_bump_gen!, gaussian_gen!, rect_gen!,
    surfgen_func_selector!,
    save_solver_config, load_solver_config, make_solver_config,


    ## Solver
    α, α0,
    M_ker, N_ker,
    M_invariant!, N_invariant!,
    solve!, solve_MDRC!



end # module RayleighSolver
