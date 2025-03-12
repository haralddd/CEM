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
    return _pre(n) / factorial(n) * (q * (q - p) * (-a0)^(n-1) + (-a0)^(n+1))
end
function _M22n(p, q, a, kperp, kpara, n)
    return _pre(n) / factorial(n) * ((q * (p - q)) * a^(n-1) / kperp - a^(n+1) / kpara)
end
function _N1n(a0, n)
    return I_pre(a0, n)
end
function _N2n(p, k, a0, n)
    return _pre(n) / factorial(n) * (k * (p - k) * a0^(n-1) - a0^(n+1))
end

function A(u::Uniaxial)
    return √((u.mu_para * u.eps_para) / (u.mu_perp * u.eps_perp))
end

function alpha_p(q::Float64, A::ComplexF64, μεpe::ComplexF64)::ComplexF64
    return A * √(μεpe - q^2)
end
function alpha_s(q::Float64, μεpa::ComplexF64)::ComplexF64
    return √(μεpa - q^2)
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

    PM11 = @view PMn[1:half, 1:half, :]
    PM12 = @view PMn[1:half, half+1:end, :]
    PM21 = @view PMn[half+1:end, 1:half, :]
    PM22 = @view PMn[half+1:end, half+1:end, :]

    PN1 = @view PNn[1:half, 1:half, :]
    PN2 = @view PNn[1:half, half+1:end, :]

    SM11 = @view SMn[1:half, 1:half, :]
    SM12 = @view SMn[1:half, half+1:end, :]
    SM21 = @view SMn[half+1:end, 1:half, :]
    SM22 = @view SMn[half+1:end, half+1:end, :]

    SN1 = @view SNn[1:half, :, :]
    SN2 = @view SNn[half+1:end, :, :]

    # --- nu = p ---
    # n = 0
    for i in axes(PM11, 1)
        p = ps[i]
        a0 = alpha0(p)
        ap = alpha_p(p, A, μεpe)
        PM11[i, i, 1] = -1
        PM12[i, i, 1] = 1
        PM21[i, i, 1] = -a0
        PM22[i, i, 1] = -ap/εpa
    end
    # n > 0
    for n in axes(PM11, 3)[2:end]
        for qidx in axes(PM11, 2)
            for pidx in axes(PM11, 1)
                p = ps[pidx]
                q = qs[qidx]
                a0 = alpha0(q)
                ap = alpha_p(q, A, μεpe)
                PM11[pidx, qidx, n] = _M11n(a0, n)
                PM12[pidx, qidx, n] = _M12n(ap, n)
                PM21[pidx, qidx, n] = _M21n(p, q, a0, n)
                PM22[pidx, qidx, n] = _M22n(p, q, ap, εpe, εpa, n)
            end
        end
    end

    for kidx in axes(PN1, 2)
        for qidx


    for pidx in axes(N1, 1)
        for kidx in axes(N1, 2)
            p = ps[pidx]
            k = ks[kidx]
            N1[pidx, kidx] = _N1n(a0, k)
            N2[pidx, kidx] = _N2n(p, k, a0)
        end
    end

    # --- nu = s ---
    # n = 0 for nu=s
    for i in axes(PM11, 1)
        p = ps[i]
        a0 = alpha0(p)
        as = alpha_s(p, μεpa)
        SM11[i, i, 1] = -1
        SM12[i, i, 1] = 1
        SM21[i, i, 1] = -a0
        SM22[i, i, 1] = -as/μpa
    end

    return nothing
end

function solve_single!(alloc::Preallocated, data::SolverData{Parameters{_S,Vacuum,Uniaxial}})::Nothing where {_S}
    
    params = data.params
    pre = data.precomputed

    FFT = params.FFT

    sFys_pqidxs = params.sFys_pqidxs
    kis = params.kis

    ys = alloc.ys
    Fys = alloc.Fys
    sFys = alloc.sFys

    PMn = pre.PMpqn
    PNn = pre.PNpkn
    PM = alloc.PMpq
    PN = alloc.PNpk

    SMn = pre.SMpqn
    SNn = pre.SNpkn
    SM = alloc.SMpq
    SN = alloc.SNpk

    PM .= 0.0
    PN .= 0.0

    SM .= 0.0
    SN .= 0.0

    @inbounds for n in reverse(axes(PMn, 3))
        for i in eachindex(Fys)
            Fys[i] = ys[i]^(n - 1)
        end

        FFT * Fys
        fftshift!(sFys, Fys)

        # Process blocks contiguously
        half_rows = size(PM, 1) ÷ 2
        
        # First process all upper blocks
        for j in 1:2:size(PM, 2)
            lj = 1 + j ÷ 2
            for i in 1:half_rows
                idx = sFys_pqidxs[i, lj]

                PM[i, j] += PMn[i, j, n] * sFys[idx]
                PM[i, j+1] += PMn[i, j+1, n] * sFys[idx]
                SM[i, j] += SMn[i, j, n] * sFys[idx]
                SM[i, j+1] += SMn[i, j+1, n] * sFys[idx]
            end
        end

        # Then process all lower blocks
        for j in 1:2:size(PM, 2)
            lj = 1 + j ÷ 2
            for i in 1:half_rows
                idx = sFys_pqidxs[i, lj]

                PM[i + half_rows, j] += PMn[i + half_rows, j, n] * sFys[idx]
                PM[i + half_rows, j+1] += PMn[i + half_rows, j+1, n] * sFys[idx]
                SM[i + half_rows, j] += SMn[i + half_rows, j, n] * sFys[idx]
                SM[i + half_rows, j+1] += SMn[i + half_rows, j+1, n] * sFys[idx]
            end
        end

        # Process PN and SN in blocks
        for j in axes(PN, 2)
            kj = kis[j]
            # Upper block
            for i in 1:half_rows
                idx = sFys_pqidxs[i, kj]
                PN[i, j] += PNn[i, j, n] * sFys[idx]
                SN[i, j] += SNn[i, j, n] * sFys[idx]
            end
            # Lower block
            for i in 1:half_rows
                idx = sFys_pqidxs[i, kj]
                PN[i + half_rows, j] += PNn[i + half_rows, j, n] * sFys[idx]
                SN[i + half_rows, j] += SNn[i + half_rows, j, n] * sFys[idx]
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