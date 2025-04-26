abstract type RandomSurface end

"Flat surface for model verification"
struct FlatSurface <: RandomSurface end

"Gaussian surface, one of the most common random surface types"
@kwdef struct GaussianSurface <: RandomSurface
    d::Float64
    a::Float64
end

"Single bump surface, for verification without simplifying to flat (close to Fresnel simplification)"
@kwdef struct SingleBumpSurface <: RandomSurface
    d::Float64
    a::Float64
end

"""
    Rectangular West O'Donnel surface
    Has very special wave number filtering properties based on cutoffs: lower, km (k⁻), and upper, kp (k⁺).
"""
@kwdef struct RectangularSurface <: RandomSurface
    d::Float64
    km::Float64
    kp::Float64
end

function pretty(surf::RandomSurface)
    io = IOBuffer()
    println(io, "$(typeof(surf))(")
    for field in fieldnames(surf)
        println(io, "\t$field=$(getfield(surf, field)),")
    end
    println(io, ")")
    return String(take!(io))
end

"""
    scale(surf::T, scale::Float64)::T where {T<:RandomSurface}

Return a new surface with the length [L] dimensions in `surf` scaled by the factor `scale`.
"""
function scale(::T, ::Float64)::T where {T<:RandomSurface}
    @error "scale not implemented for $(T)"
    return T()
end

scale(::FlatSurface, ::Float64)::FlatSurface = FlatSurface()
scale(surf::GaussianSurface, scale::Float64)::GaussianSurface = GaussianSurface(surf.d * scale, surf.a * scale)
scale(surf::SingleBumpSurface, scale::Float64)::SingleBumpSurface = SingleBumpSurface(surf.d * scale, surf.a * scale)
scale(surf::RectangularSurface, scale::Float64)::RectangularSurface = RectangularSurface(surf.d * scale, surf.km, surf.kp)