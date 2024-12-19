"""
Solves the full set of integral equations under the Rayleigh hypothesis.
I.e. solves for both R and T, so double the equations of RRE which only uses R
"""

function I(gamma, n)
    return _pre(n) * gamma^n / factorial(n)
end

function _M11n(a0, n)
    return -I(-a0, n)
end
function _M12n(a, n)
    return I(a, n)
end
function _M21n(p, q, a0, n)
    return (-1)^n * _pre(n) / factorial(n) * (q * (q - p) * a0^(n-1) + a0^(n+1))
end
function _M22n(p, q, a, kperp, kpara, n)
    return _pre(n)/(kperp * kpara * factorial(n))*(kpara * (q * (p - q)) * a^(n-1) + kperp * a^(n+1))
end
function _N1n(a0, n) # TODO: Check correct signs and such
    return I(a0, n)
end
function _N2n(p, k, a0, n) # TODO: Check correct signs and such
    return (k * (k - p) / a0 + a0) * I(a0, n)
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

function alpha(q::Float64, A::ComplexF64, μεpa::ComplexF64)::ComplexF64
    return A * (μεpa - q^2)
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

    # n = 0
    for j in axes(PMn, 2)[begin:2:end] # strided indices
        lj = 1 + j ÷ 2 # Linear idx for q
        
        for i in axes(PMn, 1)[begin:2:end] # strided indices
            li = 1 + i ÷ 2 # Linear idx for p

            p = ps[li]
            q = qs[lj]
            a = alpha(q, A, μεpa)
            a0 = alpha0(q)

            PMn[i, j, 1] = _M110(p, q)
            PMn[i, j+1, 1] = _M120(p, q)
            PMn[i+1, j, 1] = _M210(p, q, a0)
            PMn[i+1, j+1, 1] = _M220(p, q, a, εpa)

            SMn[i, j, 1] = _M110(p, q)
            SMn[i, j+1, 1] = _M120(p, q)
            SMn[i+1, j, 1] = _M210(p, q, a0)
            SMn[i+1, j+1, 1] = _M220(p, q, a, μpa)
        end
    end
    for j in axes(PNn, 2)
        for i in axes(PNn, 1)[begin:2:end]
            li = 1 + i ÷ 2
            p = ps[li]
            k = ks[j]
            a0 = alpha0(k)

            PNn[i, j, 1] = _N10(p, k)
            PNn[i+1, j, 1] = _N20(p, k, a0)

            SNn[i, j, 1] = _N10(p, k)
            SNn[i+1, j, 1] = _N20(p, k, a0)
        end
    end

    # n > 0
    for n in axes(PMn, 3)[2:end]
        for j in axes(PMn, 2)[begin:2:end]
            lj = 1 + j ÷ 2
            for i in axes(PMn, 1)[begin:2:end]
                li = 1 + i ÷ 2

                p = ps[li]
                q = qs[lj]
                a = alpha(q, A, μεpa)
                a0 = alpha0(q)

                PMn[i, j, n] = _M11n(a0, n - 1)
                PMn[i, j+1, n] = _M12n(a, n - 1)
                PMn[i+1, j, n] = _M21n(p, q, a0, n - 1)
                PMn[i+1, j+1, n] = _M22n(p, q, a, εpa, εpe, n - 1)

                SMn[i, j, n] = _M11n(a0, n - 1)
                SMn[i, j+1, n] = _M12n(a, n - 1)
                SMn[i+1, j, n] = _M21n(p, q, a0, n - 1)
                SMn[i+1, j+1, n] = _M22n(p, q, a, μpa, μpe, n - 1)
            end
        end
        for j in axes(PNn, 2)
            for i in axes(PNn, 1)[begin:2:end]
                li = 1 + i ÷ 2
                p = ps[li]
                k = ks[j]
                a0 = alpha0(k)

                PNn[i, j, n] = _N1n(a0, n - 1)
                PNn[i+1, j, n] = _N2n(p, k, a0, n - 1)

                SNn[i, j, n] = _N1n(a0, n - 1)
                SNn[i+1, j, n] = _N2n(p, k, a0, n - 1)
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

    @inbounds for n in reverse(axes(PMn, 3)) # Reverse because prefactors vanish at higher powers of ´n´
        for i in eachindex(Fys)
            Fys[i] = ys[i]^(n - 1)
        end

        FFT * Fys # In place FFT
        fftshift!(sFys, Fys)

        for j in axes(PM, 2)[begin:2:end]
            lj = 1 + j ÷ 2
            for i in axes(PM, 1)[begin:2:end]
                li = 1 + i ÷ 2
                idx = sFys_pqidxs[li, lj]

                PM[i, j] += PMn[i, j, n] * sFys[idx]
                PM[i+1, j] += PMn[i+1, j, n] * sFys[idx]
                PM[i, j+1] += PMn[i, j+1, n] * sFys[idx]
                PM[i+1, j+1] += PMn[i+1, j+1, n] * sFys[idx]

                SM[i, j] += SMn[i, j, n] * sFys[idx]
                SM[i+1, j] += SMn[i+1, j, n] * sFys[idx]
                SM[i, j+1] += SMn[i, j+1, n] * sFys[idx]
                SM[i+1, j+1] += SMn[i+1, j+1, n] * sFys[idx]
            end
        end
        for (j, kj) in enumerate(kis)
            for i in axes(PN, 1)[begin:2:end]
                li = 1 + i ÷ 2
                idx = sFys_pqidxs[li, kj]

                PN[i, j] += PNn[i, j, n] * sFys[idx]
                PN[i+1, j] += PNn[i+1, j, n] * sFys[idx]

                SN[i, j] += SNn[i, j, n] * sFys[idx]
                SN[i+1, j] += SNn[i+1, j, n] * sFys[idx]
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

"Specialization of calc_mdrc from reduced. Must split R and T of the resulting vectors."
function calc_mdrc(data::SolverData{Parameters{_S,Vacuum,Uniaxial}})::DataMDRC where {_S}
    params = data.params
    qs = params.qs
    ks = params.ks

    mask = qs .>= -1.0 .&& qs .<= 1.0
    θss = asind.(qs[mask])
    θis = data.params.θs

    @assert all(ks .≈ sind.(θis))
    Rp = data.P_res.R[begin:2:end] # Filter out R from X = (R,T)
    R2p = data.P_res.R²[begin:2:end]

    Rs = data.S_res.R[begin:2:end]
    R2s = data.S_res.R²[begin:2:end]

    coh_p, inc_p = get_coh_inc(Rp[mask, :], R2p[mask, :], qs[mask], ks)
    coh_s, inc_s = get_coh_inc(Rs[mask, :], R2s[mask, :], qs[mask], ks)

    return DataMDRC(coh_p, inc_p, coh_s, inc_s, θis, θss)
end

function calc_mdtc(data::SolverData{Parameters{_S,Vacuum,Uniaxial}}) where {_S}
end