const prefactors = [-1.0im, -1.0 + 0.0im, 1.0im, 1.0 + 0.0im] # Lookup table for imaginary prefactor
_pre(n) = prefactors[mod1(n, 4)]

"""
    function M_invariant!(M::Array{ComplexF64, 3}, parameters::Parameters{SurfT<:RandomSurface, Above<:Material, Below<:Material})

Calculates the surface invariant part of the Mpq matrix. 
Invariant under surface change and incident angle θ0/k
"""
function M_invariant!(::Array{ComplexF64,3}, params::Parameters{_S,_A,_B}, nu::Symbol)::Nothing where {_S<:RandomSurface,_A<:Material,_B<:Material}
    @error "Not implemented for $(typeof(params))"
end

"""
    function N_invariant!(N::Array{ComplexF64, 3}, parameters::Parameters{SurfT<:RandomSurface, Above<:Material, Below<:Material})

Calculates the surface invariant part of the Npk matrix. 
Invariant under surface change
"""
function N_invariant!(::Array{ComplexF64,3}, ::Parameters{_S,_A,_B}, nu::Symbol)::Nothing where {_S<:RandomSurface,_A<:Material,_B<:Material}
    @error "Not implemented for $(typeof(params))"
end

function _M_isotropic_ker(p::Float64, q::Float64, params::Parameters, kappa::ComplexF64, n::Int)::ComplexF64
    a0 = alpha(q, params.above)
    a = alpha(p, params.below)
    da = a - a0
    return _pre(n) / factorial(n) * (
        (p + kappa * q) * (p - q) * da^(n - 1) +
        (a + kappa * a0) * da^n)
end

function M_invariant!(Mpqn::Array{ComplexF64,3}, params::Parameters{_S,Vacuum,Isotropic}, nu::Symbol)::Nothing where {_S}

    ps = params.ps
    qs = params.qs
    kappa = nu == :p ? params.below.eps : params.below.mu

    # n = 0
    Mpqn[:, :, 1] .= 0.0
    @inbounds for i in axes(Mpqn, 1)
        p = ps[i]
        a0 = alpha0(p)
        a = alpha(p, params.below)
        Mpqn[i, i, 1] = a + kappa * a0
    end

    # n > 0
    @inbounds for n in axes(Mpqn, 3)[2:end]
        for j in axes(Mpqn, 2)
            for i in axes(Mpqn, 1)
                p = ps[i]
                q = qs[j]
                Mpqn[i, j, n] = _M_isotropic_ker(p, q, params, kappa, n - 1)
            end
        end
    end

    return nothing
end




function _N_isotropic_ker(p::Float64, k::Float64, params::Parameters, kappa::ComplexF64, n::Int)::ComplexF64
    a = alpha(p, params.below)
    a0 = alpha(k, params.above)
    da = a + a0
    return -_pre(n) / factorial(n) * (
        (p + kappa * k) * (p - k) * da^(n - 1) +
        (a - kappa * a0) * da^n)
end
function N_invariant!(Npkn::Array{ComplexF64,3}, params::Parameters{_S,Vacuum,Isotropic}, nu::Symbol)::Nothing where {_S}

    ps = params.ps
    ks = params.ks
    kappa = nu == :p ? params.below.eps : params.below.mu

    # n = 0
    # Npkn[:, :, 1] .= 0.0
    # @inbounds for i in axes(Npkn, 2)
    #     k = ks[i]
    #     a0 = alpha0(k)
    #     a = alpha(k, params.below)
    #     Npkn[i, i, 1] = -(a - kappa * a0)
    # end

    # n > 0
    @inbounds for n in axes(Npkn, 3)[1:end]
        for j in axes(Npkn, 2)
            for i in axes(Npkn, 1)
                p = ps[i]
                k = ks[j]
                Npkn[i, j, n] = _N_isotropic_ker(p, k, params, kappa, n - 1)
            end
        end
    end
    return nothing
end

function solve_single_reduced!(alloc::Preallocated, pre::Precomputed, data::SolverData{T})::Nothing where {T}
    @error "Not implemented for $(typeof(data))"
end



"""
    function solve_single_reduced!(alloc::Preallocated, pre::Precomputed, data::SolverData{Parameters{_S,Vacuum,Isotropic}})::Nothing where {_S}

Calculates the preallocated surface integral.
Matrix Mpqn is preallocated and contains the invariant parts of the Mpq matrix.
Matrix Npkn is preallocated and contains the invariant parts of the Npk vector.
Stores the resulting reflection coefficients in sp.Npk
Overwrites sp.Mpq with the LU factorization in the linear solution process.

# Arguments:
- `data`: [`SolverData`](@ref) - Contains the parameters, preallocated steps, output
"""
function solve_single_reduced!(alloc::Preallocated, pre::Precomputed, data::SolverData{Parameters{_S,Vacuum,Isotropic}})::Nothing where {_S}

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
