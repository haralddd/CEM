
"""
    function solve_single!(pre::Preallocated, data::SolverData)::Nothing

Calculates the preallocated surface integral.
Matrix Mpqn is preallocated and contains the invariant parts of the Mpq matrix.
Matrix Npkn is preallocated and contains the invariant parts of the Npk vector.
Stores the resulting reflection coefficients in sp.Npk
Overwrites sp.Mpq with the LU factorization in the linear solution process.

# Arguments:
- `data`: [`SolverData`](@ref) - Contains the parameters, preallocated steps, output
"""
function solve_single!(alloc::Preallocated, data::SolverData{Parameters{_S,Vacuum,Isotropic}})::Nothing where {_S}

    params = data.params
    pre = data.precomputed
    
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
    
    for n in reverse(axes(PMpqn, 3)) # Reverse because prefactors vanish at higher powers of ´n´
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

    return nothing
end
