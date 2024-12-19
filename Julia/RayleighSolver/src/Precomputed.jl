struct Precomputed
    PMpqn::Array{ComplexF64,3}
    PNpkn::Array{ComplexF64,3}

    SMpqn::Array{ComplexF64,3}
    SNpkn::Array{ComplexF64,3}

    function Precomputed(Nq::Int, Nk::Int, Ni::Int)::Precomputed
        PMpqn = Array{ComplexF64, 3}(undef, Nq, Nq, Ni)
        PNpkn = Array{ComplexF64, 3}(undef, Nq, Nk, Ni)
        
        SMpqn = similar(PMpqn)
        SNpkn = similar(PNpkn)
        return new(PMpqn, PNpkn, SMpqn, SNpkn)
    end
    function Precomputed(params::Parameters{_S,_A,_B})::Precomputed where {_S,_A,_B}
        return Precomputed(length(params.qs), length(params.ks), params.Ni+1)
    end
    function Precomputed(params::Parameters{_S,Vacuum,Uniaxial})::Precomputed where {_S}
        return Precomputed(2*length(params.qs), length(params.ks), params.Ni+1)
    end
end

"Various uniaxial crystal prefactors"
struct UniaxialParams
    kmpa::ComplexF64
    kmpe::ComplexF64
    A2::ComplexF64
    B::ComplexF64
    Cpa::ComplexF64
    Cpe::ComplexF64
end

function _A2(cr::Uniaxial)::ComplexF64
    return (cr.mu_para * cr.eps_para) / (cr.mu_perp * cr.eps_perp)
end
function UniaxialParams(below::Uniaxial, above::Vacuum, nu::Symbol)::UniParams

    if nu == :p
        # k = κ     p/m = ±     pa/pe = ∥/⟂
        kmpa = below.eps_para
        kmpe = below.eps_perp
        Cpa = kmpa
        Cpe = kmpe
    else
        # k = κ     p/m = ±     pa/pe = ∥/⟂
        kmpa = below.mu_para
        kmpe = below.mu_perp
        Cpa = kmpa
        Cpe = kmpe
    end
    A2 = _A2(below)
    B = kmpa / kmpe
    return UniParams(kmpa, kmpe, A2, B, Cpa, Cpe)
end

function alpha2(q::Float64, A2::ComplexF64, με::ComplexF64)::ComplexF64
    return A2 * (με - q^2)
end
function alpha2(q::Float64, material::Uniaxial)::ComplexF64
    A2 = _A2(material)
    με = material.mu_para * material.eps_para
    return alpha2(q, A2, με)
end
function alpha_tilde(q, p, params::UniaxialParams)::ComplexF64
    B = params.B
    A2 = params.A2
    με = params.με
    return sqrt(B * (q^2 - p^2) + alpha2(q, A2, με))
end

"""
    alpha(q::Float64, material<:Material)::ComplexF64
Calculates the perpendicular wave number, alpha ≡ q⟂, based on the parallel wave number `q` and `material` properties.

# Arguments:
- `q`: Wavenumber - Dimensionless, scaled by ``\\frac{\\omega}{c}``
- `material`: `Material` containing information about the material
"""
function alpha(q::Float64, material::T) where {T<:Material} end
function alpha(q::Float64, material::Vacuum)::ComplexF64
    return sqrt(complex(1.0 - q^2))
end
function alpha(q::Float64, material::Isotropic)::ComplexF64
    με = material.eps * material.mu
    return sqrt(με - q^2)
end
function alpha(q::Float64, material::Uniaxial)::ComplexF64
    return sqrt(alpha2(q, material))
end
function alpha0(q)
    alpha(q, Vacuum())
end

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
    @inbounds for j in axes(Mpqn, 2)
        for i in axes(Mpqn, 1)
            p = ps[i]
            q = qs[j]
            a0 = alpha0(q)
            a = alpha(p, params.below)
            Mpqn[i, j, 1] = (p ≈ q) ? a + kappa * a0 : 0.0
        end
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


function _M_uniaxial_ker(p::Float64, q::Float64, params::Parameters, upa::UniaxialParams, n::Int)::ComplexF64
    kmpa = upa.kmpa
    kmpe = upa.kmpe
    Cpa = upa.Cpa
    Cpe = upa.Cpe

    at = alpha_tilde(q, p, upa)
    a0 = alpha(q, params.above)
    da = at - a0
    return _pre(n) / factorial(n) * (
        kmpa * (p + Cpe * q) * (p - q) * da^(n - 1) +
        kmpe * (at + Cpa * a0) * da^n)
end

function M_invariant!(Mpqn::Array{ComplexF64,3}, params::Parameters{_S,Vacuum,Uniaxial}, nu::Symbol)::Nothing where {_S}

    ps = params.ps
    qs = params.qs
    upa = UniaxialParams(params)

    @inbounds for n in axes(Mpqn, 3)
        for j in axes(Mpqn, 2)
            for i in axes(Mpqn, 1)
                p = ps[i]
                q = qs[j]
                Mpqn[i, j, n] = _M_uni_ker(p, q, n - 1, params, upa)
            end
        end
    end
    return nothing
end


"""
    function N_invariant!(N::Array{ComplexF64, 3}, parameters::Parameters{SurfT<:RandomSurface, Above<:Material, Below<:Material})

Calculates the surface invariant part of the Npk matrix. 
Invariant under surface change
"""
function N_invariant!(::Array{ComplexF64,3}, ::Parameters{_S,_A,_B}, nu::Symbol)::Nothing where {_S<:RandomSurface,_A<:Material,_B<:Material}
    @error "Not implemented for $(typeof(params))"
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
    @inbounds for j in axes(Npkn, 2)
        for i in axes(Npkn, 1)
            p = ps[i]
            k = ks[j]
            a0 = alpha0(k)
            a = alpha(p, params.below)
            Npkn[i, j, 1] = (p ≈ k) ? a - kappa * a0 : 0.0
        end
    end

    # n > 0
    @inbounds for n in axes(Npkn, 3)[2:end]
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

function _N_uniaxial_ker(p::Float64, q::Float64, params::Parameters, upa::UniaxialParams, n::Int)::ComplexF64
    kmpa = upa.kmpa
    kmpe = upa.kmpe
    Cpa = upa.Cpa
    Cpe = upa.Cpe

    at = alpha_tilde(q, p, upa)
    a0 = alpha(q, params.above)
    da = at + a0
    return -_pre(n) / factorial(n) * (
        kmpa * (p + Cpe * q) * (p - q) * da^(n - 1) +
        kmpe * (at - Cpa * a0) * da^n)
end

function N_invariant!(N::Array{ComplexF64,3}, params::Parameters{_S,Vacuum,Uniaxial}, nu::Symbol)::Nothing where {_S}

    ps = params.ps
    qs = params.qs
    kis = params.kis

    upa = UniaxialParams(params)


    @inbounds for n in axes(N, 3), (j, kj) in enumerate(kis), i in axes(N, 1)
        p = ps[i]
        k = qs[kj]
        N[i, j, n] = _N_uni_ker(p, k, params, upa, n - 1)
    end
    return nothing
end

function precompute!(pre::Precomputed, params::Parameters{S,A,B})::Nothing where {S,A,B}
    M_invariant!(pre.PMpqn, params, :p)
    N_invariant!(pre.PNpkn, params, :p)

    M_invariant!(pre.SMpqn, params, :s)
    N_invariant!(pre.SNpkn, params, :s)
    return nothing
end

function _validate_single(A)
    if !all(isfinite.(A))
        open("trace.dat", "w+") do io
            write(io, A)
        end
        error("Validation of precomputed data failed, matrices not finite")
    end
end

function validate(pre::Precomputed)::Nothing
    _validate_single(pre.PMpqn)
    _validate_single(pre.PNpkn)
    _validate_single(pre.SMpqn)
    _validate_single(pre.SNpkn)
    return nothing
end