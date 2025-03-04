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
    return _pre(n) / factorial(n) * (q * (q - p) * a0^(n-1) + a0^(n+1))
end
function _M22n(p, q, a, kperp, kpara, n)
    return _pre(n)/factorial(n) * 1/(kperp * kpara) * (kpara * (q * (q - p)) * a^(n-1) + kperp * a^(n+1))
end
function _N1n(a0, n) # TODO: Check correct signs and such
    return -I(-a0, n)
end
function _N2n(p, k, a0, n) # TODO: Check correct signs and such
    return _pre(n) / factorial(n) * (k * (k - p) * a0^(n-1) + a0^(n+1))
end

# n=0 specializations, I(n=0) = δ(p-q) from the Fourier integral, since ζ^0 = 1
δ(x) = (x ≈ 0.0) ? 1.0 : 0.0  # Allow for machine epsilon
function _M110(p, q)
    return -δ(p - q)
end
function _M120(p, q)
    return δ(p - q)
end
function _M210(p, q, a0)
    return δ(p - q) * a0
end
function _M220(p, q, a, kpara)
    return δ(p - q) * a / kpara
end
function _N10(p, k)
    return δ(p - k)
end
function _N20(p, k, a0)
    return δ(p - k) * a0
end

function A(u::Uniaxial)
    return √((u.mu_para * u.eps_para) / (u.mu_perp * u.eps_perp))
end

function alpha(q::Float64, A::ComplexF64, μεpa::ComplexF64)::ComplexF64
    return A * √(μεpa - q^2)
end

function alpha(q::Float64, uni::Uniaxial)::ComplexF64
    _A = A(uni)
    μεpa = uni.mu_para * uni.eps_para
    return alpha(q, _A, μεpa)
end

function precompute!(pre::Precomputed, params::Parameters{_S,Vacuum,Uniaxial})::Nothing where {_S}
    εpe = params.below.eps_perp
    εpa = params.below.eps_para
    μpe = params.below.mu_perp
    μpa = params.below.mu_para
    A = √((μpa * εpa) / (μpe * εpe))
    μεpa = μpa * εpa

    qs = params.qs
    ps = params.ps
    ks = params.ks

    PMn = pre.PMpqn
    PNn = pre.PNpkn
    SMn = pre.SMpqn
    SNn = pre.SNpkn

    # Separate system of equations into 
    # n = 0
    for j in 1:2:size(PMn, 2) # Process column pairs
        lj = 1 + j ÷ 2 # Linear idx for q
        
        # First block: all M11 and M12 elements
        for i in 1:2:size(PMn, 1)
            li = 1 + i ÷ 2 # Linear idx for p
            p = ps[li]
            q = qs[lj]
            a = alpha(q, A, μεpa)
            a0 = alpha0(q)

            PMn[li, j, 1] = _M110(p, q)
            PMn[li, j+1, 1] = _M120(p, q)
            SMn[li, j, 1] = _M110(p, q)
            SMn[li, j+1, 1] = _M120(p, q)
        end

        # Second block: all M21 and M22 elements
        for i in 1:2:size(PMn, 1)
            li = 1 + i ÷ 2
            p = ps[li]
            q = qs[lj]
            a = alpha(q, A, μεpa)
            a0 = alpha0(q)

            PMn[li + size(PMn,1)÷2, j, 1] = _M210(p, q, a0)
            PMn[li + size(PMn,1)÷2, j+1, 1] = _M220(p, q, a, εpa)
            SMn[li + size(SMn,1)÷2, j, 1] = _M210(p, q, a0)
            SMn[li + size(SMn,1)÷2, j+1, 1] = _M220(p, q, a, μpa)
        end
    end
    for j in axes(PNn, 2)
        # First block: all N1 elements
        for i in 1:2:size(PNn, 1)
            li = 1 + i ÷ 2
            p = ps[li]
            k = ks[j]
            a0 = alpha0(k)

            PNn[li, j, 1] = _N10(p, k)
            SNn[li, j, 1] = _N10(p, k)
        end

        # Second block: all N2 elements
        for i in 1:2:size(PNn, 1)
            li = 1 + i ÷ 2
            p = ps[li]
            k = ks[j]
            a0 = alpha0(k)

            PNn[li + size(PNn,1)÷2, j, 1] = _N20(p, k, a0)
            SNn[li + size(SNn,1)÷2, j, 1] = _N20(p, k, a0)
        end
    end

    # n > 0
    for n in axes(PMn, 3)[2:end]
        for j in 1:2:size(PMn, 2)
            lj = 1 + j ÷ 2

            # First block: M11 and M12 elements
            for i in 1:2:size(PMn, 1)
                li = 1 + i ÷ 2
                p = ps[li]
                q = qs[lj]
                a = alpha(q, A, μεpa)
                a0 = alpha0(q)

                PMn[li, j, n] = _M11n(a0, n - 1)
                PMn[li, j+1, n] = _M12n(a, n - 1)
                SMn[li, j, n] = _M11n(a0, n - 1)
                SMn[li, j+1, n] = _M12n(a, n - 1)
            end

            # Second block: M21 and M22 elements
            for i in 1:2:size(PMn, 1)
                li = 1 + i ÷ 2
                p = ps[li]
                q = qs[lj]
                a = alpha(q, A, μεpa)
                a0 = alpha0(q)

                PMn[li + size(PMn,1)÷2, j, n] = _M21n(p, q, a0, n - 1)
                PMn[li + size(PMn,1)÷2, j+1, n] = _M22n(p, q, a, εpa, εpe, n - 1)
                SMn[li + size(SMn,1)÷2, j, n] = _M21n(p, q, a0, n - 1)
                SMn[li + size(SMn,1)÷2, j+1, n] = _M22n(p, q, a, μpa, μpe, n - 1)
            end
        end

        for j in axes(PNn, 2)
            # First block: N1 elements
            for i in 1:2:size(PNn, 1)
                li = 1 + i ÷ 2
                p = ps[li]
                k = ks[j]
                a0 = alpha0(k)

                PNn[li, j, n] = _N1n(a0, n - 1)
                SNn[li, j, n] = _N1n(a0, n - 1)
            end

            # Second block: N2 elements
            for i in 1:2:size(PNn, 1)
                li = 1 + i ÷ 2
                p = ps[li]
                k = ks[j]
                a0 = alpha0(k)

                PNn[li + size(PNn,1)÷2, j, n] = _N2n(p, k, a0, n - 1)
                SNn[li + size(SNn,1)÷2, j, n] = _N2n(p, k, a0, n - 1)
            end
        end
    end
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