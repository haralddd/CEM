function energy_ratio(R2, T2, ks, qs, params::Parameters{_S,Vacuum,Uniaxial}, ν=:p) where {_S}
    ret = zeros(size(R2, 2))
    μpa = params.below.mu_para
    μpe = params.below.mu_perp
    εpa = params.below.eps_para
    εpe = params.below.eps_perp
    μεpa = μpa * εpa
    Aval = A(params.below)
    μεpe = μpe * εpe

    κpa = ν == :p ? εpa : μpa

    for k_idx in eachindex(ks)  # For each incident angle
        k = ks[k_idx]
        a0k = alpha0(k)
        for (q_idx, q) in enumerate(qs)
            a0q = alpha0(q)
            aq = ν == :p ? alpha_p(q, Aval, μεpe) : alpha_s(q, μεpa)

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

function energy_conservation(P_res::Results, S_res::Results, params::Parameters{_S,Vacuum,Uniaxial}) where {_S}

    ks = params.ks
    qs = params.qs

    filt = -1.0 .<= qs .<= 1.0
    R2p = get_R²(P_res)[filt, :]
    T2p = get_T²(P_res)[filt, :]
    R2s = get_R²(S_res)[filt, :]
    T2s = get_T²(S_res)[filt, :]

    P = energy_ratio(R2p, T2p, ks, qs[filt], params, :p)
    S = energy_ratio(R2s, T2s, ks, qs[filt], params, :s)

    return P, S
end

function energy_conservation(PR2::Matrix{Float64}, SR2::Matrix{Float64}, params::Parameters{_S,Vacuum,Isotropic}) where {_S}
    ks = params.ks
    qs = params.qs

    filt = qs .>= -1.0 .&& qs .<= 1.0
    R2p = PR2[filt, :]
    R2s = SR2[filt, :]

    P = energy_ratio(R2p, ks, qs[filt])
    S = energy_ratio(R2s, ks, qs[filt])

    return P, S
end

function energy_conservation(P_res::Results, S_res::Results, params::Parameters{_S,Vacuum,Isotropic}) where {_S}
    energy_conservation(P_res.R², S_res.R², params)
end