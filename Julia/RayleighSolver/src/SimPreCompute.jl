"""

    struct SimPreCompute

Contains the precomputed parts of the system of equations. 
This amounts to the factors which are invariant under surface change
"""
struct SimPreCompute
    Mpqn::Array{ComplexF64,3}
    Npkn::Array{ComplexF64,3}
end

"""
    function M_invariant!(M::Array{ComplexF64, 3}, parameters::SimParams{SurfT<:RandomSurface, PolarizationT<:Polarization, Above<:Material, Below<:Material})

Calculates the surface invariant part of the Mpq matrix. 
Invariant under surface change and incident angle θ0/k
"""
function M_invariant!(::Array{ComplexF64,3}, params::SimParams{_S,_P,_A,_B})::Nothing where {_S<:RandomSurface,_P<:Polarization,_A<:Material,_B<:Material}
    @error "Not implemented for $(typeof(params))"
end

"""
    function N_invariant!(N::Array{ComplexF64, 3}, parameters::SimParams{SurfT<:RandomSurface, PolarizationT<:Polarization, Above<:Material, Below<:Material})

Calculates the surface invariant part of the Npk matrix. 
Invariant under surface change
"""
function N_invariant!(::Array{ComplexF64,3}, ::SimParams{_S,_P,_A,_B})::Nothing where {_S<:RandomSurface,_P<:Polarization,_A<:Material,_B<:Material}
    @error "Not implemented for $(typeof(params))"
end

const prefactors = [-1.0im, -1.0 + 0.0im, 1.0im, 1.0 + 0.0im] # Lookup table
_pre(n) = prefactors[mod1(n, 4)]

function M_invariant!(M::Array{ComplexF64,3}, params::SimParams{_S,Nu,Vacuum,Isotropic})::Nothing where {_S,Nu<:Polarization}
    
    ps = params.ps
    qs = params.qs
    above = params.above
    below = params.below
    kappa = typeof(Nu) == PolarizationP ? below.eps : below.mu

    function _M_iso_ker(p::Float64, q::Float64, n::Int)::ComplexF64
        a0 = alpha(q, above)
        a = alpha(p, below)
        da = a - a0
        return _pre(n) / factorial(n) * (
            (p + kappa * q) * (p - q) * da^(-1) +
            (a + kappa * a0))*da^n
    end

    @inbounds for n in axes(M, 3), j in axes(M, 2), i in axes(M, 1)
        p = ps[i]
        q = qs[j]
        M[i, j, n] = _M_iso_ker(p, q, n - 1)
    end
    return nothing
end

function N_invariant!(N::Array{ComplexF64,3}, params::SimParams{_S,Nu,Vacuum,Isotropic})::Nothing where {_S,Nu<:Polarization}

    ps = params.ps
    qs = params.qs
    kis = params.kis
    above = params.above
    below = params.below
    
    kappa = typeof(Nu) == PolarizationP ? below.eps : below.mu

    function _N_iso_ker(p::Float64, k::Float64, n::Int)::ComplexF64
        a0 = alpha(k, above)
        a = alpha(p, below)
        da = a + a0
        return _pre(n) / factorial(n) * (
            (p + kappa * k) * (p - k) * da^(n - 1) +
            (a - kappa * a0) * da^n
        )
    end

    @inbounds for n in axes(N, 3), (j, kj) in enumerate(kis), i in axes(N, 1)
        p = ps[i]
        k = qs[kj]
        N[i, j, n] = _N_iso_ker(p, k, n - 1)
    end
    return nothing
end

"Bnu ≡ μ∥/μ⟂ for nu=p and ε∥/ε⟂ for nu=s"
function ptilde(q::Float64, p::Float64, kmt::ComplexF64)::ComplexF64
    return √(
        (1 - kmt) * q^2 + kmt * p^2
    )
end

function M_invariant!(M::Array{ComplexF64,3}, params::SimParams{_S,Nu,UniaxialCrystal,UniaxialCrystal})::Nothing where {_S,Nu<:Polarization}

    ps = params.ps
    qs = params.qs
    above::UniaxialCrystal = params.above
    below::UniaxialCrystal = params.below
    if typeof(Nu) == PolarizationP
        # k = κ     p/m = ±     pa/pe = ∥/⟂
        kmpa = below.eps_para
        kmpe = below.eps_perp
        kppa = above.eps_para
        kppe = above.eps_perp
        kmt = below.mu_para / below.mu_perp # kappa minus tilde, used in ptilde
        Bpe = kmpe / kppe
        Bpa = kmpa / kppa
    else
        kmpa = below.mu_para
        kmpe = below.mu_perp
        kppa = above.mu_para
        kppe = above.mu_perp
        kmt = below.eps_para / below.eps_perp
        Bpe = kmpe / kppe
        Bpa = kmpa / kppa
    end

    function _M_uni_ker(p::Float64, q::Float64, n::Int)::ComplexF64
        pt = ptilde(q,p,kmt)
        a = alpha(p, below)
        a0 = alpha(q, above)
        da = a - a0
        return _pre(n) / factorial(n) * kmpa * (
            (pt + Bpe * q) * (pt - q) * da^(-1) +
            kmpe * (a + Bpa * a0)
        ) * da^n
    end

    @inbounds for n in axes(M, 3), j in axes(M, 2), i in axes(M, 1)
        p = ps[i]
        q = qs[j]
        M[i, j, n] = _M_uni_ker(p, q, n - 1)
    end
    return nothing
end

function N_invariant!(N::Array{ComplexF64,3}, params::SimParams{_S,Nu,UniaxialCrystal,UniaxialCrystal})::Nothing where {_S,Nu<:Polarization}

    ps = params.ps
    qs = params.qs
    kis = params.kis
    above::UniaxialCrystal = params.above
    below::UniaxialCrystal = params.below
    if typeof(Nu) == PolarizationP
        # k = κ     p/m = ±     pa/pe = ∥/⟂
        kmpa = below.eps_para
        kmpe = below.eps_perp
        kppa = above.eps_para
        kppe = above.eps_perp
        kmt = below.mu_para / below.mu_perp # kappa minus tilde, used in ptilde
        Bpe = kmpe / kppe
        Bpa = kmpa / kppa
    else
        kmpa = below.mu_para
        kmpe = below.mu_perp
        kppa = above.mu_para
        kppe = above.mu_perp
        kmt = below.eps_para / below.eps_perp
        Bpe = kmpe / kppe
        Bpa = kmpa / kppa
    end

    function _N_uni_ker(p::Float64, k::Float64, n::Int)::ComplexF64
        pt = ptilde(k,p,kmt)
        a0 = alpha(k, above)
        a = alpha(p, below)
        da = a + a0
        return _pre(n) / factorial(n) * kmpa * (
            (pt + Bpe * k) * (pt - k) * da^(- 1) +
            kmpe * (a - Bpa * a0)
        ) * da^n
    end

    @inbounds for n in axes(N, 3), (j, kj) in enumerate(kis), i in axes(N, 1)
        p = ps[i]
        k = qs[kj]
        N[i, j, n] = _N_uni_ker(p, k, n - 1)
    end
    return nothing
end

function M_invariant(rp::SimParams{_S,Nu,A,B})::Array{ComplexF64,3} where {_S,Nu,A,B}
    Mpqn = Array{ComplexF64,3}(undef, length(rp.ps), length(rp.qs), rp.Ni + 1)
    M_invariant!(Mpqn, rp)
    return Mpqn
end
function N_invariant(rp::SimParams{_S,Nu,A,B})::Array{ComplexF64,3} where {_S,Nu,A,B}
    Npkn = Array{ComplexF64,3}(undef, length(rp.ps), length(rp.kis), rp.Ni + 1)
    N_invariant!(Npkn, rp)
    return Npkn
end

function SimPreCompute(rp::SimParams)::SimPreCompute
    return SimPreCompute(M_invariant(rp), N_invariant(rp))
end

function validate(pc::SimPreCompute)::Nothing
    @assert all(isfinite.(pc.Mpqn))
    @assert all(isfinite.(pc.Npkn))
    return nothing
end