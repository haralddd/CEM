#! Must be included after simulation_parameters.jl and simulation_prealloc.jl

""" 
    correlation(k::Float64, surf<:RandomSurface)::Float64

Correlation function of the RandomSurface in Fourier space"""
function correlation(::Float64, ::T)::Float64 where {T<:RandomSurface} end

function correlation(k::Float64, surf::GaussianSurface)::Float64
    a = surf.a
    return √π * a * exp(-(a * k / 2.0)^2)
end

function correlation(k::Float64, surf::RectangularSurface)::Float64
    kp = surf.kp
    km = surf.km
    H(x) = (x > 0.0) ? 1.0 : 0.0

    return π / (kp - km) * (H(kp - k) * H(k - km) + H(kp + k) * H(-k - km))
end

"""
    generate_surface!(sp::SimPrealloc, parameters::SimParams{SurfT,_MA,_MB})

Generate random surface of the given type and overwrites sp.ys in-place
Reverts to generic Fourier filtering function, which requires correlation(k, spa.surf::T) to be implemented
""" 
function generate_surface!(sp::SimPrealloc, spa::SimParams{SurfT,_MA,_MB})::Nothing where {SurfT <: RandomSurface, _MA, _MB}
    d = spa.surf.d
    xks = spa.xks
    rng = spa.rng
    FFT = spa.FFT
    IFFT = spa.IFFT
    A = 1/sqrt(spa.dx)

    ys = sp.ys
    Z = sp.Z

    # Generate random phases
    randn!(rng, ys)
    @inbounds for i in eachindex(Z)
        Z[i] = complex(ys[i] * d)
    end
    FFT * Z # In place FFT
    @inbounds for i in eachindex(ys)
        ys[i] = sqrt(correlation(xks[i], spa.surf))*A
        Z[i] = Z[i] * ys[i] # Filter the numbers in Fourier space
    end
    IFFT * Z # in place Inverse FFT
    @inbounds for i in eachindex(ys)
        ys[i] = real(Z[i])
    end

    return nothing
end
function generate_surface!(sp::SimPrealloc, ::SimParams{FlatSurface,_MA,_MB})::Nothing where {_MA,_MB}
    sp.ys .= 0.0
    return nothing
end
function generate_surface!(sp::SimPrealloc, spa::SimParams{SingleBumpSurface,_MA,_MB})::Nothing where {_MA,_MB}
    surf = spa.surf
    δ = surf.d
    a = surf.a

    sp.ys .= δ * exp.((-0.5 * spa.xs / a) .^ 2)
    return nothing
end