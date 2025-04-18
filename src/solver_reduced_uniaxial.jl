
function _M_uniaxial_ker(p::Float64, q::Float64, kappa_para::ComplexF64, kappa_perp::ComplexF64, alpha_func::Function, n::Int)::ComplexF64
    a = alpha_func(p)
    a0 = alpha(q, params.above)
    da = a - a0

    return _pre(n) / factorial(n) * (
        kappa_para * (p + kappa_perp * q) * (p - q) * da^(n - 1) +
        kappa_perp * (a + kappa_para * a0) * da^n
    )
end

function _N_uniaxial_ker(p::Float64, k::Float64, kappa_para::ComplexF64, kappa_perp::ComplexF64, alpha_func::Function, n::Int)::ComplexF64
    a = alpha_func(p)
    a0 = alpha0(k)
    da = a + a0

    return -_pre(n) / factorial(n) * (
        kappa_para * (p + kappa_perp * k) * (p - k) * da^(n - 1) +
        kappa_perp * (a - kappa_para * a0) * da^n
    )
end

"""
    M_invariant_uniaxial!(Mpqn::Array{ComplexF64,3}, params::Parameters{_S,Vacuum,Uniaxial}, nu::Symbol)::Nothing where {_S}

Implementation of M_invariant_uniaxial! for Vacuum-Uniaxial interface.

# Arguments:
- `Mpqn`: 3D array to store the precomputed matrix elements
- `params`: [`Parameters`](@ref) containing simulation parameters
- `nu`: Polarization, either `:p` or `:s`
"""
function M_invariant!(Mpqn::Array{ComplexF64,3}, params::Parameters{_S,Vacuum,Uniaxial}, nu::Symbol)::Nothing where {_S}
    ps = params.ps
    qs = params.qs
    below = params.below

    @assert params.mu_para == params.mu_perp "Uniaxial material must have isotropic permeability to use Reduced Uniaxial solver"

    # Define parameters based on polarization
    kappa_para = nu == :p ? below.eps_para : below.mu_para
    kappa_perp = nu == :p ? below.eps_perp : below.mu_perp
    alpha_func = nu == :p ? (p -> alpha_p(p, below)) : (p -> alpha_s(p, below))

    # n = 0
    Mpqn[:, :, 1] .= 0.0
    @inbounds for i in axes(Mpqn, 1)
        p = ps[i]
        a0 = alpha0(p)
        a = alpha_func(p)

        # Diagonal term for n=0
        Mpqn[i, i, 1] = kappa_perp * (a + kappa_para * a0)
    end

    # n > 0
    @inbounds for n in axes(Mpqn, 3)[2:end]
        for j in axes(Mpqn, 2)
            for i in axes(Mpqn, 1)
                p = ps[i]
                q = qs[j]
                Mpqn[i, j, n] = _M_uniaxial_ker(p, q, kappa_para, kappa_perp, alpha_func, n - 1)
            end
        end
    end

    return nothing
end


"""
    N_invariant!(Npkn::Array{ComplexF64,3}, params::Parameters{_S,Vacuum,Uniaxial}, nu::Symbol)::Nothing where {_S}

Implementation of N_invariant! for Vacuum-Uniaxial interface.

# Arguments:
- `Npkn`: 3D array to store the precomputed matrix elements
- `params`: [`Parameters`](@ref) containing simulation parameters
- `nu`: Polarization, either `:p` or `:s`
"""
function N_invariant!(Npkn::Array{ComplexF64,3}, params::Parameters{_S,Vacuum,Uniaxial}, nu::Symbol)::Nothing where {_S}
    ps = params.ps
    ks = params.ks
    below = params.below

    # Define parameters based on polarization
    kappa_para = nu == :p ? below.eps_para : below.mu_para
    kappa_perp = nu == :p ? below.eps_perp : below.mu_perp
    alpha_func = nu == :p ? (p -> alpha_p(p, below)) : (p -> alpha_s(p, below))

    # Process all n values
    @inbounds for n in axes(Npkn, 3)
        for j in axes(Npkn, 2)
            for i in axes(Npkn, 1)
                p = ps[i]
                k = ks[j]
                Npkn[i, j, n] = _N_uniaxial_ker(p, k, kappa_para, kappa_perp, alpha_func, n - 1)
            end
        end
    end

    return nothing
end

"""
    function solve_single!(alloc::Preallocated, pre::Precomputed, data::SolverData{Parameters{_S,Vacuum,Uniaxial}})::Nothing where {_S}

Calculates the surface integral for a given surface realization in a uniaxial material.
Uses the precomputed invariant parts of M and N matrices.

# Arguments:
- `alloc`: [`Preallocated`](@ref) structure containing working arrays
- `pre`: [`Precomputed`](@ref) structure containing precomputed values
- `data`: [`SolverData`](@ref) structure containing the parameters and results
"""
function solve_single!(alloc::Preallocated, pre::Precomputed, data::SolverData{Parameters{_S,Vacuum,Uniaxial}})::Nothing where {_S}
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

    return nothing
end
