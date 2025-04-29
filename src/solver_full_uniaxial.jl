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

"""
    M_invariant_full!(Mpqn::Array{ComplexF64,3}, params::Parameters{ST,Vacuum,BT})::Nothing where {ST,BT}

Compute the M invariant matrices for both P and S polarizations for the full solver.
"""
function M_invariant_full!(Mpqn::Array{ComplexF64,3}, params::Parameters{ST,Vacuum,BT}, nu::Symbol)::Nothing where {ST,BT}
    below = params.below
    
    if below isa Uniaxial
        εpe = below.eps_perp
        εpa = below.eps_para
        μpe = below.mu_perp
        μpa = below.mu_para
        A = get_A(below)
        kappa_perp = nu == :p ? εpe : μpe
        kappa_para = nu == :p ? εpa : μpa
        alpha_func = nu == :p ? (q -> alpha_p(q, below)) : (q -> alpha_s(q, below))
    elseif below isa Isotropic
        εpe = εpa = below.eps
        μpe = μpa = below.mu
        A = 1.0 + 0.0im
        alpha_func = (q -> alpha(q, below))
        kappa_perp = kappa_para = nu == :p ? εpe : μpe
    else
        error("Full solver only supports Vacuum-Uniaxial and Vacuum-Isotropic interfaces")
    end

    qs = params.qs
    ps = params.ps

    Mpqn .= 0.0+0.0im
    
    half = length(ps)
    @assert half == size(Mpqn, 1) ÷ 2
    
    Mpqn11 = @view Mpqn[1:half, 1:half, :]
    Mpqn12 = @view Mpqn[1:half, half+1:end, :]
    Mpqn21 = @view Mpqn[half+1:end, 1:half, :]
    Mpqn22 = @view Mpqn[half+1:end, half+1:end, :]
    
    # M-elements
    # n = 0
    for i in axes(Mpqn11, 1)
        p = ps[i]
        a0 = alpha0(p)
        a = alpha_func(p)
        Mpqn11[i, i, 1] = -1
        Mpqn12[i, i, 1] = 1
        Mpqn21[i, i, 1] = -a0
        Mpqn22[i, i, 1] = -a / kappa_para
    end
    # n > 0
    for n in axes(Mpqn11, 3)[2:end]
        for qidx in axes(Mpqn11, 2)
            for pidx in axes(Mpqn11, 1)
                p = ps[pidx]
                q = qs[qidx]
                a0 = alpha0(q)
                a = alpha_func(q)
                Mpqn11[pidx, qidx, n] = _M11n(a0, n-1)
                Mpqn12[pidx, qidx, n] = _M12n(a, n-1)
                Mpqn21[pidx, qidx, n] = _M21n(p, q, a0, n-1)
                Mpqn22[pidx, qidx, n] = _M22n(p, q, a, kappa_perp, kappa_para, n-1)
            end
        end
    end
    return nothing
end

"""
    N_invariant_full!(Npkn::Array{ComplexF64,3}, params::Parameters{ST,Vacuum,BT})::Nothing where {ST,BT}

Compute the N invariant matrices for both P and S polarizations for the full solver.
"""
function N_invariant_full!(Npkn::Array{ComplexF64,3}, params::Parameters{ST,Vacuum,BT}, nu::Symbol)::Nothing where {ST,BT}
    below = params.below

    if below isa Uniaxial
        kappa_para = nu == :p ? below.eps_para : below.mu_para
        kappa_perp = nu == :p ? below.eps_perp : below.mu_perp
    elseif below isa Isotropic
        kappa_perp = kappa_para = nu == :p ? below.eps : below.mu
    else
        error("Uniaxial solver only supports Vacuum-Uniaxial and Vacuum-Isotropic interfaces")
    end
    
    ps = params.ps
    ks = params.ks

    Npkn .= 0.0+0.0im

    half = length(ps)
    @assert half == size(Npkn, 1) ÷ 2
    
    Npkn1 = @view Npkn[1:half, :, :]
    Npkn2 = @view Npkn[half+1:end, :, :]
    
    # N-elements
    for n in axes(Npkn1, 3)
        for kidx in axes(Npkn1, 2)
            for pidx in axes(Npkn1, 1)
                p = ps[pidx]
                k = ks[kidx]
                a0 = alpha0(k)
                Npkn1[pidx, kidx, n] = _N1n(a0, n-1)
                Npkn2[pidx, kidx, n] = _N2n(p, k, a0, n-1)
            end
        end
    end
    
    return nothing
end

function solve_single_full!(alloc::Preallocated, pre::Precomputed, data::SolverData)::Nothing
    @error "Not implemented for $(typeof(data))"
end

function solve_single_full!(alloc::Preallocated, pre::Precomputed, data::SolverData{Parameters{ST,Vacuum,BT}})::Nothing where {ST,BT}

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