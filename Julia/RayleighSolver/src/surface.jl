"Surface parameter superclass for different SurfType surface types"
abstract type SurfaceParams end
struct FlatSurfaceParams <: SurfaceParams end
struct GaussianSurfaceParams <: SurfaceParams
    d::Float64
    a::Float64
end
struct SingleBumpSurfaceParams <: SurfaceParams
    d::Float64
    a::Float64
end
struct RectSurfaceParams <: SurfaceParams
    d::Float64
    km::Float64
    kp::Float64
end


"""
    scale(surf::T, scale::Float64)::T where {T<:SurfaceParams}

Scales the physical parameters in `surf` by the scaling factor `scale`.
"""
function scale(surface::T, scale::Float64)::T where {T<:SurfaceParams}
    @error "scale not implemented for $(T)"
    return T()
end

scale(::FlatSurfaceParams, ::Float64)::FlatSurfaceParams = FlatSurfaceParams()
scale(surf::GaussianSurfaceParams, scale::Float64)::GaussianSurfaceParams = GaussianSurfaceParams(surf.d * scale, surf.a * scale)
scale(surf::SingleBumpSurfaceParams, scale::Float64)::SingleBumpSurfaceParams = SingleBumpSurfaceParams(surf.d * scale, surf.a * scale)
scale(surf::RectSurfaceParams, scale::Float64)::RectSurfaceParams = RectSurfaceParams(surf.d * scale, surf.km, surf.kp)

