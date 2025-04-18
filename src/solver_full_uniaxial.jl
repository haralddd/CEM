"""
Solves the full set of integral equations under the Rayleigh hypothesis.
I.e. solves for both R and T, so double the equations of RRE which only uses R
"""

function I_pre(gamma, n)
    return _pre(n) * gamma^n / factorial(n)
end

function _M11n(a0, n)
    return -I_pre(-a0, n)
end
function _M12n(a, n)
    return I_pre(a, n)
end
function _M21n(p, q, a0, n)
    return _pre(n) / factorial(n) * (q * (q - p) * (-a0)^(n - 1) + (-a0)^(n + 1))
end
function _M22n(p, q, a, kperp, kpara, n)
    return _pre(n) / factorial(n) * ((q * (p - q)) * a^(n - 1) / kperp - a^(n + 1) / kpara)
end
function _N1n(a0, n)
    return I_pre(a0, n)
end
function _N2n(p, k, a0, n)
    return _pre(n) / factorial(n) * (k * (p - k) * a0^(n - 1) - a0^(n + 1))
end



function precompute!(pre::Precomputed, params::Parameters{_S,Vacuum,Uniaxial})::Nothing where {_S}
    εpe = params.below.eps_perp
    εpa = params.below.eps_para
    μpe = params.below.mu_perp
    μpa = params.below.mu_para
    A = √((μpa * εpa) / (μpe * εpe))
    μεpa = μpa * εpa
    μεpe = μpe * εpe

    qs = params.qs
    ps = params.ps
    ks = params.ks

    PMn = pre.PMpqn
    PNn = pre.PNpkn
    SMn = pre.SMpqn
    SNn = pre.SNpkn

    PMn .= 0.0
    PNn .= 0.0
    SMn .= 0.0
    SNn .= 0.0

    half = length(ps)
    @assert half == size(PMn, 1) ÷ 2 == size(SMn, 1) ÷ 2 == size(PNn, 1) ÷ 2 == size(SNn, 1) ÷ 2

    PMn11 = @view PMn[1:half, 1:half, :]
    PMn12 = @view PMn[1:half, half+1:end, :]
    PMn21 = @view PMn[half+1:end, 1:half, :]
    PMn22 = @view PMn[half+1:end, half+1:end, :]

    PNn1 = @view PNn[1:half, :, :]
    PNn2 = @view PNn[half+1:end, :, :]

    SMn11 = @view SMn[1:half, 1:half, :]
    SMn12 = @view SMn[1:half, half+1:end, :]
    SMn21 = @view SMn[half+1:end, 1:half, :]
    SMn22 = @view SMn[half+1:end, half+1:end, :]

    SNn1 = @view SNn[1:half, :, :]
    SNn2 = @view SNn[half+1:end, :, :]

    # --- nu = p ---
    # M-elements
    # n = 0
    for i in axes(PMn11, 1)
        p = ps[i]
        a0 = alpha0(p)
        ap = alpha_p(p, A, μεpa)
        PMn11[i, i, 1] = -1
        PMn12[i, i, 1] = 1
        PMn21[i, i, 1] = -a0
        PMn22[i, i, 1] = -ap / εpa
    end
    # n > 0
    for n in axes(PMn11, 3)[2:end]
        for qidx in axes(PMn11, 2)
            for pidx in axes(PMn11, 1)
                p = ps[pidx]
                q = qs[qidx]
                a0 = alpha0(q)
                ap = alpha_p(q, A, μεpa)
                PMn11[pidx, qidx, n] = _M11n(a0, n-1)
                PMn12[pidx, qidx, n] = _M12n(ap, n-1)
                PMn21[pidx, qidx, n] = _M21n(p, q, a0, n-1)
                PMn22[pidx, qidx, n] = _M22n(p, q, ap, εpe, εpa, n-1)
            end
        end
    end
    # N-elements
    for n in axes(PMn11, 3)
        for kidx in axes(PNn1, 2)
            for pidx in axes(PNn1, 1)
                p = ps[pidx]
                k = ks[kidx]
                a0 = alpha0(k)
                PNn1[pidx, kidx, n] = _N1n(a0, n-1)
                PNn2[pidx, kidx, n] = _N2n(p, k, a0, n-1)
            end
        end
    end

    # --- nu = s ---
    # M-elements
    # n = 0
    for i in axes(SMn11, 1)
        p = ps[i]
        a0 = alpha0(p)
        as = alpha_s(p, μεpa)
        SMn11[i, i, 1] = -1
        SMn12[i, i, 1] = 1
        SMn21[i, i, 1] = -a0
        SMn22[i, i, 1] = -as / μpa
    end
    # n > 0
    for n in axes(SMn11, 3)[2:end]
        for qidx in axes(SMn11, 2)
            for pidx in axes(SMn11, 1)
                p = ps[pidx]
                q = qs[qidx]
                a0 = alpha0(q)
                as = alpha_s(q, μεpa)
                SMn11[pidx, qidx, n] = _M11n(a0, n-1)
                SMn12[pidx, qidx, n] = _M12n(as, n-1)
                SMn21[pidx, qidx, n] = _M21n(p, q, a0, n-1)
                SMn22[pidx, qidx, n] = _M22n(p, q, as, μpe, μpa, n-1)
            end
        end
    end
    # N-elements
    for n in axes(SNn1, 3)
        for kidx in axes(SNn1, 2)
            for pidx in axes(SNn1, 1)
                p = ps[pidx]
                k = ks[kidx]
                a0 = alpha0(k)
                SNn1[pidx, kidx, n] = _N1n(a0, n-1)
                SNn2[pidx, kidx, n] = _N2n(p, k, a0, n-1)
            end
        end
    end

    return nothing
end

function solve_single!(alloc::Preallocated, pre::Precomputed, data::SolverData{Parameters{_S,Vacuum,Uniaxial}})::Nothing where {_S}

    params = data.params

    FFT = params.FFT

    sFys_pqidxs = params.sFys_pqidxs
    kis = params.kis

    ys = alloc.ys
    Fys = alloc.Fys
    sFys = alloc.sFys

    half = length(params.ps)

    # Partial matrices
    PMn = pre.PMpqn
    PNn = pre.PNpkn
    SMn = pre.SMpqn
    SNn = pre.SNpkn

    PMn11 = @view PMn[1:half, 1:half, :]
    PMn12 = @view PMn[1:half, half+1:end, :]
    PMn21 = @view PMn[half+1:end, 1:half, :]
    PMn22 = @view PMn[half+1:end, half+1:end, :]

    PNn1 = @view PNn[1:half, :, :]
    PNn2 = @view PNn[half+1:end, :, :]

    SMn11 = @view SMn[1:half, 1:half, :]
    SMn12 = @view SMn[1:half, half+1:end, :]
    SMn21 = @view SMn[half+1:end, 1:half, :]
    SMn22 = @view SMn[half+1:end, half+1:end, :]

    SNn1 = @view SNn[1:half, :, :]
    SNn2 = @view SNn[half+1:end, :, :]

    # Solution matrices
    PM = alloc.PMpq
    PN = alloc.PNpk
    SM = alloc.SMpq
    SN = alloc.SNpk

    PM .= 0.0
    PN .= 0.0
    SM .= 0.0
    SN .= 0.0

    PM11 = @view PM[1:half, 1:half]
    PM12 = @view PM[1:half, half+1:end]
    PM21 = @view PM[half+1:end, 1:half]
    PM22 = @view PM[half+1:end, half+1:end]

    PN1 = @view PN[1:half, :]
    PN2 = @view PN[half+1:end, :]

    SM11 = @view SM[1:half, 1:half]
    SM12 = @view SM[1:half, half+1:end]
    SM21 = @view SM[half+1:end, 1:half]
    SM22 = @view SM[half+1:end, half+1:end]

    SN1 = @view SN[1:half, :]
    SN2 = @view SN[half+1:end, :]

    for n in reverse(axes(PMn11, 3))
        for i in eachindex(Fys)
            Fys[i] = ys[i]^(n - 1)
        end

        FFT * Fys
        fftshift!(sFys, Fys)

        for j in axes(PM11, 2)
            for i in axes(PM11, 1)
                idx = sFys_pqidxs[i, j]
                PM11[i, j] += PMn11[i, j, n] * sFys[idx]
                PM12[i, j] += PMn12[i, j, n] * sFys[idx]
                PM21[i, j] += PMn21[i, j, n] * sFys[idx]
                PM22[i, j] += PMn22[i, j, n] * sFys[idx]

                SM11[i, j] += SMn11[i, j, n] * sFys[idx]
                SM12[i, j] += SMn12[i, j, n] * sFys[idx]
                SM21[i, j] += SMn21[i, j, n] * sFys[idx]
                SM22[i, j] += SMn22[i, j, n] * sFys[idx]
            end
        end
        for (j, kj) in enumerate(kis)
            for i in axes(PN1, 1)
                idx = sFys_pqidxs[i, kj]
                PN1[i, j] += PNn1[i, j, n] * sFys[idx]
                PN2[i, j] += PNn2[i, j, n] * sFys[idx]

                SN1[i, j] += SNn1[i, j, n] * sFys[idx]
                SN2[i, j] += SNn2[i, j, n] * sFys[idx]
            end
        end
    end

    A = lu!(PM)
    @inbounds for j in axes(PN, 2)
        b = @view PN[:, j]
        ldiv!(A, b)
    end

    A = lu!(SM)
    @inbounds for j in axes(SN, 2)
        b = @view SN[:, j]
        ldiv!(A, b)
    end
end