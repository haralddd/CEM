

function A_ker(a,n)
    return (-1.0im * a)^n / factorial(n)
end

function B1_ker(a0,n)
    return (-1.0im * a0)^n / factorial(n)
end

function B2_ker(a0,n)
    return (1.0im * a0)^n / factorial(n)
end

function A_invariant_hybrid!(A::Matrix{ComplexF64}, params::Parameters{ST,Vacuum,BT}, nu::Symbol)::Nothing where {ST,BT}
    below = params.below
    if below isa Uniaxial
        alpha_func = nu == :p ? (q -> alpha_p(q, below)) : (q -> alpha_s(q, below))
    elseif below isa Isotropic
        alpha_func = (q -> alpha(q, below))
    else
        error("Hybrid solver only supports Vacuum-Uniaxial and Vacuum-Isotropic interfaces")
    end

    qs = params.qs
    
    # A-elements
    for n in axes(A, 2)
        for qidx in axes(A, 1)
            q = qs[qidx]
            a = alpha_func(q)
            A[qidx, n] = A_ker(a, n-1)
        end
    end
    
    return nothing
end

function B1_invariant_hybrid!(B1::Matrix{ComplexF64}, params::Parameters{ST,Vacuum,BT})::Nothing where {ST,BT}
    ks = params.ks
    # B1-elements
    for n in axes(B1, 2)
        for kidx in axes(B1, 1)
            k = ks[kidx]
            a0 = alpha0(k)
            B1[kidx, n] = B1_ker(a0, n-1)
        end
    end
    
    return nothing
end

function B2_invariant_hybrid!(B2::Matrix{ComplexF64}, params::Parameters{ST,Vacuum,BT})::Nothing where {ST,BT}
    qs = params.qs
    # B2-elements
    for n in axes(B2, 2)
        for qidx in axes(B2, 1)
            q = qs[qidx]
            a0 = alpha0(q)
            B2[qidx, n] = B2_ker(a0, n-1)
        end
    end
    
    return nothing
end

function calc_T!(alloc::Preallocated, precomputed::Precomputed, params::Parameters)
    qs = params.qs
    ps = params.ps
    ks = params.ks
    kis = params.kis

    ys = alloc.ys
    Fys = alloc.Fys
    sFys = alloc.sFys
    FFT = params.FFT
    sFys_pqidxs = params.sFys_pqidxs
    

    PApq = alloc.hybrid.PApq
    Pbpk = alloc.hybrid.Pbpk
    SApq = alloc.hybrid.SApq
    Sbpk = alloc.hybrid.Sbpk

    PAqn = precomputed.hybrid.PAqn
    PB1kn = precomputed.hybrid.PB1kn
    PB2qn = precomputed.hybrid.PB2qn
    SAqn = precomputed.hybrid.SAqn
    SB1kn = precomputed.hybrid.SB1kn
    SB2qn = precomputed.hybrid.SB2qn

    for n in reverse(axes(PAqn, 2))
        for i in eachindex(Fys)
            Fys[i] = ys[i]^(n - 1)
        end
        FFT * Fys
        fftshift!(sFys, Fys)

        for (kidx, kqidx) in enumerate(kis)
            for pidx in eachindex(ps)
                fourier_idx = sFys_pqidxs[pidx, kqidx]

                Pbpk[pidx, kidx] += PB1kn[kidx, n] * sFys[fourier_idx]
                Sbpk[pidx, kidx] += SB1kn[kidx, n] * sFys[fourier_idx]
            end

            for qidx in eachindex(qs)
                for pidx in eachindex(ps)
                    fourier_idx = sFys_pqidxs[pidx, qidx]

                    PApq[pidx, qidx] += PAqn[qidx, n] * sFys[fourier_idx]
                    SApq[pidx, qidx] += SAqn[qidx, n] * sFys[fourier_idx]
                    Pbpk[pidx, kidx] += PB2qn[qidx, n] * sFys[fourier_idx]
                    Sbpk[pidx, kidx] += SB2qn[qidx, n] * sFys[fourier_idx]
                end
            end
        end
    end

    # Solve the linear systems
    A = lu!(PApq)
    @inbounds for j in axes(Pbpk, 2)
        b = @view Pbpk[:, j]
        ldiv!(A, b)
    end

    A = lu!(SApq)
    @inbounds for j in axes(Sbpk, 2)
        b = @view Sbpk[:, j]
        ldiv!(A, b)
    end
    return nothing
end



"""
    function solve_single_hybrid!(alloc::Preallocated, pre::Precomputed, data::SolverData{Parameters{ST,Vacuum,BT}})::Nothing where {ST,BT}

Calculates the surface integral for a given surface realization in a uniaxial material using the hybrid method.
Uses the precomputed invariant parts of M, N, A, B1 and B2 matrices.

# Arguments:
- `alloc`: [`Preallocated`](@ref) structure containing working arrays
- `pre`: [`Precomputed`](@ref) structure containing precomputed values
- `data`: [`SolverData`](@ref) structure containing the parameters and results
"""
function solve_single_hybrid!(alloc::Preallocated, pre::Precomputed, data::SolverData{Parameters{ST,Vacuum,BT}})::Nothing where {ST,BT}
    params = data.params

    FFT = params.FFT

    sFys_pqidxs = params.sFys_pqidxs
    kis = params.kis

    ys = alloc.ys
    Fys = alloc.Fys
    sFys = alloc.sFys

    PMpqn = pre.PMpqn
    PNpkn = pre.PNpkn
    PMpq = alloc.PMpq
    PNpk = alloc.PNpk

    SMpqn = pre.SMpqn
    SNpkn = pre.SNpkn
    SMpq = alloc.SMpq
    SNpk = alloc.SNpk

    PMpq .= 0.0
    PNpk .= 0.0

    SMpq .= 0.0
    SNpk .= 0.0

    for n in reverse(axes(PMpqn, 3)) # Reverse because prefactors vanish at higher powers of 'n'
        for i in eachindex(Fys)
            Fys[i] = ys[i]^(n - 1)
        end

        FFT * Fys # In place FFT
        fftshift!(sFys, Fys)

        for j in axes(PMpq, 2)
            for i in axes(PMpq, 1)
                idx = sFys_pqidxs[i, j]
                PMpq[i, j] += PMpqn[i, j, n] * sFys[idx]
                SMpq[i, j] += SMpqn[i, j, n] * sFys[idx]
            end
        end

        for (j, kj) in enumerate(kis)
            for i in axes(PNpk, 1)
                idx = sFys_pqidxs[i, kj]
                PNpk[i, j] += PNpkn[i, j, n] * sFys[idx]
                SNpk[i, j] += SNpkn[i, j, n] * sFys[idx]
            end
        end
    end

    # Solve the linear systems
    A = lu!(PMpq)
    @inbounds for j in axes(PNpk, 2)
        b = @view PNpk[:, j]
        ldiv!(A, b)
    end

    A = lu!(SMpq)
    @inbounds for j in axes(SNpk, 2)
        b = @view SNpk[:, j]
        ldiv!(A, b)
    end

    # Calculate T in place from the values of R
    calc_T!(alloc, pre, params)

    return nothing
end