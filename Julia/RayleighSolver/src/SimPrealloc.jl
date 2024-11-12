struct SystemPreAlloc
    Mpqn::Array{ComplexF64,3}
    Npkn::Array{ComplexF64,3}
    Mpq::Matrix{ComplexF64}
    Npk::Matrix{ComplexF64}

    function SystemPreAlloc(Nq::Int, Nk::Int, Ni::Int)::SystemPreAlloc
        Mpqn = zeros(ComplexF64, Nq, Nq, Ni)
        Npkn = zeros(ComplexF64, Nq, Nk, Ni)
        Mpq = zeros(ComplexF64, Nq, Nq)
        Npk = zeros(ComplexF64, Nq, Nk)
        return new(Mpqn, Npkn, Mpq, Npk)
    end
end


"""
    SimPrealloc{N<:Integer}(Mpq::Matrix{ComplexF64}, Npk::Matrix{ComplexF64}, Fys::Vector{ComplexF64}, sFys::Vector{ComplexF64}, ys::Vector{Float64})

SimPrealloc{N<:Integer} is a struct which contains all
preallocated matrices which are used during
the calculation of the surface integral, i.e. the array values are mutated in place
Note that all members are uninitialized, and must be initialized after making the struct

# Fields

- `Mpq::Matrix{ComplexF64}`: Matrix of the Mpq coefficients (A)
- `Npk::Matrix{ComplexF64}`: Vector (for all k) of Vectors of the Npk coefficients (b)
- `Fys::Vector{ComplexF64}`: Fourier transform of surface heights, prealloc
- `sFys::Vector{ComplexF64}`: Shifted Fourier transform of surface heights, prealloc
- `Z::Vector{ComplexF64}`: Preallocated computation step in surface generation
- `ys::Vector{Float64}`: Surface heights
"""
struct SimPrealloc
    p_data::SystemPreAlloc
    s_data::SystemPreAlloc
    Fys::Vector{ComplexF64}
    sFys::Vector{ComplexF64}
    ys::Vector{Float64}

    function SimPrealloc(Nq::Int, Nk::Int, Ni::Int)
        p_data = SystemPreAlloc(Nq, Nk, Ni)
        s_data = SystemPreAlloc(Nq, Nk, Ni)

        ys = Vector{Float64}(undef, 2Nq)
        Fys = similar(ys, ComplexF64)
        sFys = similar(Fys)

        new(p_data, s_data, Fys, sFys, ys)
    end
    function SimPrealloc(spa::SimParams)::SimPrealloc
        SimPrealloc(spa.Nq, length(spa.ks), spa.Ni)
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

function _A2(cr::UniaxialCrystal)::ComplexF64
    return (cr.mu_para*cr.eps_para)/(cr.mu_perp*cr.eps_perp)
end
function UniaxialParams(below::UniaxialCrystal, above::Vacuum, nu::Symbol)::UniParams
    
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
function alpha2(q::Float64, material::UniaxialCrystal)::ComplexF64
    A2 = _A2(material)
    με = material.mu_para*material.eps_para
    return alpha2(q, A2, με)
end
function alpha_tilde(q, p, params::UniaxialParams)::ComplexF64
    B = params.B
    A2 = params.A2
    με = params.με
    return sqrt(B*(q^2 - p^2) + alpha2(q, A2, με))
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
    με = material.eps*material.mu
    return sqrt(με - q^2)
end
function alpha(q::Float64, material::UniaxialCrystal)::ComplexF64
    return sqrt(alpha2(q, material))
end
function alpha0(q) alpha(q, Vacuum()) end

const prefactors = [-1.0im, -1.0 + 0.0im, 1.0im, 1.0 + 0.0im] # Lookup table for imaginary prefactor
_pre(n) = prefactors[mod1(n, 4)]


"""
    function M_invariant!(M::Array{ComplexF64, 3}, parameters::SimParams{SurfT<:RandomSurface, Above<:Material, Below<:Material})

Calculates the surface invariant part of the Mpq matrix. 
Invariant under surface change and incident angle θ0/k
"""
function M_invariant!(::Array{ComplexF64,3}, spa::SimParams{_S,_A,_B}, nu::Symbol)::Nothing where {_S<:RandomSurface,_A<:Material,_B<:Material}
    @error "Not implemented for $(typeof(spa))"
end

function _M_isotropic_ker(p::Float64, q::Float64, spa::SimParams, kappa::ComplexF64, n::Int)::ComplexF64
    a0 = alpha(q, spa.above)
    a = alpha(p, spa.below)
    da = a - a0
    return _pre(n) / factorial(n) * (
        (p + kappa * q) * (p - q) * da^(n-1) +
        (a + kappa * a0)*da^n)
end

function M_invariant!(Mpqn::Array{ComplexF64,3}, spa::SimParams{_S,Vacuum,Isotropic}, nu::Symbol)::Nothing where {_S}
    
    ps = spa.ps
    qs = spa.qs
    kappa = nu==:p ? spa.below.eps : spa.below.mu

    Mpqn .= 0.0
    @inbounds for i in axes(Mpqn, 1)
        p = ps[i]
        a = alpha(p, spa.below)
        a0 = alpha0(p)
        Mpqn[i, i, 1] = a + kappa*a0
    end

    @inbounds for n in axes(Mpqn, 3)[2:end], j in axes(Mpqn, 2), i in axes(Mpqn, 1)
        p = ps[i]
        q = qs[j]
        Mpqn[i, j, n] = _M_isotropic_ker(p, q, spa, kappa, n-1)
    end
    
    return nothing
end


function _M_uniaxial_ker(p::Float64, q::Float64, spa::SimParams, upa::UniaxialParams, n::Int)::ComplexF64
    kmpa = upa.kmpa
    kmpe = upa.kmpe
    Cpa = upa.Cpa
    Cpe = upa.Cpe

    at = alpha_tilde(q, p, upa)
    a0 = alpha(q, spa.above)
    da = at - a0
    return _pre(n) / factorial(n) * (
        kmpa * (p + Cpe * q) * (p - q) * da^(n-1) +
        kmpe * (at + Cpa * a0) * da^n)
end

function M_invariant!(Mpqn::Array{ComplexF64,3}, spa::SimParams{_S,Vacuum,UniaxialCrystal}, nu::Symbol)::Nothing where {_S}

    ps = spa.ps
    qs = spa.qs
    upa = UniaxialParams(spa)

    @inbounds for n in axes(Mpqn, 3), j in axes(Mpqn, 2), i in axes(Mpqn, 1)
        p = ps[i]
        q = qs[j]
        Mpqn[i, j, n] = _M_uni_ker(p, q, n - 1, spa, upa)
    end
    return nothing
end


"""
    function N_invariant!(N::Array{ComplexF64, 3}, parameters::SimParams{SurfT<:RandomSurface, Above<:Material, Below<:Material})

Calculates the surface invariant part of the Npk matrix. 
Invariant under surface change
"""
function N_invariant!(::Array{ComplexF64,3}, ::SimParams{_S,_A,_B}, nu::Symbol)::Nothing where {_S<:RandomSurface,_A<:Material,_B<:Material}
    @error "Not implemented for $(typeof(spa))"
end


function _N_isotropic_ker(p::Float64, k::Float64, spa::SimParams, kappa::ComplexF64, n::Int)::ComplexF64
    a0 = alpha(k, spa.above)
    a = alpha(p, spa.below)
    da = a + a0
    return _pre(n) / factorial(n) * (
        (p + kappa * k) * (p - k) * da^(n - 1) +
        (a - kappa * a0) * da^n)
end
function N_invariant!(Npkn::Array{ComplexF64,3}, spa::SimParams{_S,Vacuum,Isotropic}, nu::Symbol)::Nothing where {_S}

    ps = spa.ps
    ks = spa.ks
    kis = spa.kis
    
    kappa = nu==:p ? spa.below.eps : spa.below.mu
    Npkn .= 0.0
    @inbounds for (i, ki) in enumerate(kis)
        k = ks[i]

        a = alpha(k, spa.below)
        a0 = alpha0(k)
        Npkn[ki, i, 1] = a - kappa * a0
    end

    @inbounds for n in axes(Npkn, 3)[2:end], (j, k) in enumerate(ks), i in axes(Npkn, 1)
        p = ps[i]
        Npkn[i, j, n] = _N_isotropic_ker(p, k, spa, kappa, n - 1)
    end
    return nothing
end

function _N_uniaxial_ker(p::Float64, q::Float64, spa::SimParams, upa::UniaxialParams, n::Int)::ComplexF64
    kmpa = upa.kmpa
    kmpe = upa.kmpe
    Cpa = upa.Cpa
    Cpe = upa.Cpe

    at = alpha_tilde(q, p, upa)
    a0 = alpha(q, spa.above)
    da = at + a0
    return _pre(n) / factorial(n) * (
        kmpa * (p + Cpe * q) * (p - q) * da^(n-1) +
        kmpe * (at - Cpa * a0) * da^n)
end

function N_invariant!(N::Array{ComplexF64,3}, spa::SimParams{_S,Vacuum,UniaxialCrystal}, nu::Symbol)::Nothing where {_S}

    ps = spa.ps
    qs = spa.qs
    kis = spa.kis

    upa = UniaxialParams(spa)


    @inbounds for n in axes(N, 3), (j, kj) in enumerate(kis), i in axes(N, 1)
        p = ps[i]
        k = qs[kj]
        N[i, j, n] = _N_uni_ker(p, k, spa, upa, n-1)
    end
    return nothing
end

function precompute!(sp::SimPrealloc, spa::SimParams{S,A,B})::Nothing where {S,A,B}
    M_invariant!(sp.p_data.Mpqn, spa, :p)
    M_invariant!(sp.s_data.Mpqn, spa, :s)
    N_invariant!(sp.s_data.Npkn, spa, :s)
    N_invariant!(sp.p_data.Npkn, spa, :p)
    return nothing
end

function validate(sp::SimPrealloc)::Nothing
    @assert all(isfinite.(sp.p_data.Mpqn))
    @assert all(isfinite.(sp.s_data.Mpqn))
    @assert all(isfinite.(sp.p_data.Npkn))
    @assert all(isfinite.(sp.s_data.Npkn))
    return nothing
end