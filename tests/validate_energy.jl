using RayleighSolver

function scattering_condition(res::Results, ks, qs, A::ComplexF64, μεpa::ComplexF64, κpa::ComplexF64)
    filt = -1.0 .<= qs .<= 1.0
    R = get_R(res)[filt, :]
    T = get_T(res)[filt, :]
    ret = zeros(size(R, 2))

    for k_idx in eachindex(ks)  # For each incident angle
        k = ks[k_idx]
        a0k = real(alpha0(k))
        for (q_idx, q) in enumerate(qs[filt])
            a0q = real(alpha0(q))
            aq = real(alpha(q, A, μεpa))

            R_term = abs2(R[q_idx, k_idx]) * a0q / a0k
            T_term = abs2(T[q_idx, k_idx]) * aq / (real(κpa) * a0k)

            ret[k_idx] += R_term + T_term
        end
    end
    
    return ret
end


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
function validate_energy_conservation(data::SolverData{Parameters{_S,Vacuum,Uniaxial}}) where {_S}
    μpa = data.params.below.mu_para
    μpe = data.params.below.mu_perp
    εpa = data.params.below.eps_para
    εpe = data.params.below.eps_perp
    μεpa = μpa * εpa
    Aval = A(data.params.below)
    ks = data.params.ks
    qs = data.params.qs

    @info "μ∥ = $μpa"
    @info "ε∥ = $εpa"
    @info "μ⟂ = $μpe"
    @info "ε⟂ = $εpe"

    P = scattering_condition(data.P_res, ks, qs, Aval, μεpa, εpa)
    S = scattering_condition(data.S_res, ks, qs, Aval, μεpa, μpa)

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
    data_path = joinpath(@__DIR__, "..", "output")
    output_files = readdir(data_path)
    file = output_files[findfirst(x->occursin(filename, x), output_files)]
    data = load_solver_data(joinpath(data_path, file))

    cond_P, cond_S = validate_energy_conservation(data)
end
    