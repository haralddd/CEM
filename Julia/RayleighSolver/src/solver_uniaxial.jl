abstract type Anisotropy{T} end

struct Uniaxial{T} <: Anisotropy{T}
    para::T
    perp::T
end

abstract type Crystal{T} end

struct UniaxialCrystal{T} <: Crystal{T}
    eps::Uniaxial{T}
    mu::Uniaxial{T}
end

struct IsotropicCrystal{T} <: Crystal{T}
    eps::T
    mu::T
end
function A(cr::UniaxialCrystal{T})::T where T
    return √((cr.mu.para*cr.eps.para)/(cr.mu.perp*cr.eps.perp))
end

function alpha(q::Float64, cr::UniaxialCrystal{T})::T where T
    return A(cr) * sqrt(cr.mu.perp*cr.eps.perp - q^2)
end

function pt(q::Float64, p::Float64, anti_kappa::Uniaxial{T})::T where T
    b = anti_kappa.para / anti_kappa.perp
    return √(
        (1-b)*q^2 + b*p^2
    )
end