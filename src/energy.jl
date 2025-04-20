function energy_ratio(R2, T2, ks, qs, params::Parameters, ν=:p)
    ret = zeros(size(R2, 2))
    μpa = params.below.mu_para
    μpe = params.below.mu_perp
    εpa = params.below.eps_para
    εpe = params.below.eps_perp
    μεpa = μpa * εpa
    Aval = A(params.below)

    κpa = ν == :p ? εpa : μpa

    for k_idx in eachindex(ks)  # For each incident angle
        k = ks[k_idx]
        a0k = alpha0(k)
        for (q_idx, q) in enumerate(qs)
            a0q = alpha0(q)
            aq = ν == :p ? alpha_p(q, Aval, μεpa) : alpha_s(q, μεpa)

            R_term = R2[q_idx, k_idx] * a0q / a0k
            T_term = T2[q_idx, k_idx] * real(aq / κpa) / a0k

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