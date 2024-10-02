
"""
    abstract type Material

Superclass for material properties
Subclasses contains permittivity and permeability of the material
"""
abstract type Material end

struct Vacuum <: Material end

struct Isotropic <: Material
    eps::ComplexF64
    mu::ComplexF64
end

struct UniaxialCrystal <: Material
    eps_perp::ComplexF64
    eps_para::ComplexF64
    mu_perp::ComplexF64
    mu_para::ComplexF64
end

function _A(cr::UniaxialCrystal)::ComplexF64
    return √((cr.mu_para*cr.eps_para)/(cr.mu_perp*cr.eps_perp))
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
    A = _A(material)
    με = material.mu_perp*material.eps_perp
    return A * sqrt(με - q^2)
end

function alpha0(q) alpha(q, Vacuum()) end