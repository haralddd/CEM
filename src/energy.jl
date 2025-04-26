function energy_ratio(R2, T2, ks, qs, params::Parameters, ν=:p)
    ret = zeros(size(R2, 2))
    below = params.below
    if below isa Uniaxial
        μpa = below.mu_para
        μpe = below.mu_perp
        εpa = below.eps_para
        εpe = below.eps_perp
        alpha_func = ν == :p ? (q -> alpha_p(q, below)) : (q -> alpha_s(q, below))
    elseif below isa Isotropic
        μpa = μpe = below.mu
        εpa = εpe = below.eps
        alpha_func = (q -> alpha(q, below))
    end

    κpa = ν == :p ? εpa : μpa

    for k_idx in eachindex(ks)  # For each incident angle
        k = ks[k_idx]
        a0k = alpha0(k)
        for (q_idx, q) in enumerate(qs)
            a0q = alpha0(q)
            aq = alpha_func(q)

            R_term = real(R2[q_idx, k_idx] * a0q / a0k)
            T_term = real(T2[q_idx, k_idx] * aq / (κpa * a0k))

            ret[k_idx] += R_term + T_term
        end
    end
    
    return ret
end

function energy_ratio(R2, ks, qs)
    ret = zeros(size(R2, 2))

    for k_idx in eachindex(ks)  # For each incident angle
        k = ks[k_idx]
        a0k = alpha0(k)
        for (q_idx, q) in enumerate(qs)
            a0q = alpha0(q)

            ret[k_idx] += R2[q_idx, k_idx] * a0q / a0k
        end
    end
    
    return ret
end

function energy_conservation(data::SolverData)
    params = data.params
    

    ks = params.ks
    qs = params.qs

    if data.solver_type == :reduced
        filt = -1.0 .<= qs .<= 1.0
        R2p = data.Rp.A²[filt, :]
        R2s = data.Rs.A²[filt, :]

        P = energy_ratio(R2p, ks, qs[filt])
        S = energy_ratio(R2s, ks, qs[filt])
    elseif data.solver_type == :full
        R2p = data.Rp.A²
        T2p = data.Tp.A²
        R2s = data.Rs.A²
        T2s = data.Ts.A²

        P = energy_ratio(R2p, T2p, ks, qs, params, :p)
        S = energy_ratio(R2s, T2s, ks, qs, params, :s)
    end

    return P, S
end