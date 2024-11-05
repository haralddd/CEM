
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
