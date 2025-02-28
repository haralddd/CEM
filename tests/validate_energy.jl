using RayleighSolver

"""
    validate_energy_conservation(data::SolverData; rtol=1e-6)

Validates the energy conservation relation for reflection and transmission coefficients
from the solver data. The equation being validated is:
1/L₁ ∫ dq/2π [|R'(q|k)|²α０(q,ω)/α０(k,ω) + |T'(q|k)|²α(q,ω)/α０(k,ω)] = 1

Parameters:
- data: SolverData containing the solution matrices
- rtol: Relative tolerance for comparison with 1

Returns:
- Dictionary with validation results for each incident angle and polarization
"""
function validate_energy_conservation(data::SolverData)
    P, S = energy_conservation(data.P_res, data.S_res, data.params)
    ks = data.params.ks

    println("P")
    for p in P
        println(p)
    end
    println("S")
    for s in S
        println(s)
    end

    P_cond = P .≈ 1.0
    S_cond = S .≈ 1.0

    if any(.~P_cond)
        loc = findall(.~P_cond)
        @warn "P energy not conserved at k=$(ks[loc])"
    end

    if any(.~S_cond)
        loc = findall(.~S_cond)
        @warn "S energy not conserved at k=$(ks[loc])"
    end

    return P_cond, S_cond
end


if (abspath(PROGRAM_FILE) == @__FILE__) || isinteractive()
    filename = ARGS[1]
    data_path = joinpath(splitdir(@__DIR__)[1], "output")
    output_files = readdir(data_path)
    file = output_files[findfirst(x->occursin(filename, x), output_files)]
    data = load_solver_data(joinpath(data_path, file))

    cond_P, cond_S = validate_energy_conservation(data)
end
    