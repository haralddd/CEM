
"""
    abstract type Material

Superclass for material properties
Subclasses contains permittivity and permeability of the material
"""
abstract type Material end

struct Vacuum <: Material end

@kwdef struct Isotropic <: Material
    eps::ComplexF64
    mu::ComplexF64
end

@kwdef struct Uniaxial <: Material
    eps_perp::ComplexF64
    eps_para::ComplexF64
    mu_perp::ComplexF64
    mu_para::ComplexF64
end

function pretty(mat::Material)
    io = IOBuffer()
    println(io, "$(typeof(mat))(")
    for field in fieldnames(mat)
        println(io, "\t$field=$(getfield(mat, field)),")
    end
    println(io, ")")
    return String(take!(io))
end

"""
    signum(x::Float64)
Returns the sign of `x` as a float.
"""
signum(x::Float64)::Float64 = (x < 0.0) ? -1.0 : 1.0

"""
    alpha(q::Float64, material::T)::ComplexF64
Calculates the perpendicular wave number, alpha ≡ q⟂, based on the parallel wave number `q` and `material` properties.

# Arguments:
- `q`: Wavenumber - Dimensionless, scaled by ``\\frac{\\omega}{c}``
- `material`: `Material` containing information about the material
"""
function alpha(q::Float64, material::T) where {T<:Material} end
function alpha(q::Float64, material::Vacuum)::ComplexF64
    D = sqrt(complex(1.0 - q^2))
    return D*signum(imag(D))
end
function alpha(q::Float64, material::Isotropic)::ComplexF64
    με = material.eps * material.mu
    D = sqrt(με - q^2)
    # Choose the branch where the imaginary part is positive
    return D*signum(imag(D))
end
function alpha0(q)
    alpha(q, Vacuum())
end

function get_A(mat::Uniaxial)::ComplexF64
    return (mat.mu_para * mat.eps_para) / (mat.mu_perp * mat.eps_perp)
end

function alpha_p(q::Float64, A::ComplexF64, μεpa::ComplexF64)::ComplexF64
    D = sqrt(μεpa - A * q^2)
    # Choose the branch where the imaginary part is positive
    return D*signum(imag(D))
end
function alpha_s(q::Float64, μεpa::ComplexF64)::ComplexF64
    D = sqrt(μεpa - q^2)
    # Choose the branch where the imaginary part is positive
    return D*signum(imag(D))
end
function alpha_p(q::Float64, mat::Uniaxial)::ComplexF64
    _A = get_A(mat)
    μεpa = mat.mu_para*mat.eps_para
    return alpha_p(q, _A, μεpa)
end
function alpha_s(q::Float64, mat::Uniaxial)::ComplexF64
    μεpa = mat.mu_para*mat.eps_para
    return alpha_s(q, μεpa)
end