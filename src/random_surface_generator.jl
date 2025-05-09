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
    generate_surface!(sp::SimPrealloc, parameters::Parameters{SurfT,_MA,_MB})

Generate random surface of the given type and overwrites sp.ys in-place
Reverts to generic Fourier filtering function, which requires correlation(k, params.surf::T) to be implemented
"""
function generate_surface!(pre::Preallocated, params::Parameters{SurfT,_MA,_MB})::Nothing where {SurfT<:RandomSurface,_MA,_MB}
    d = params.surf.d
    xks = params.xks
    rng = params.rng
    FFT = params.FFT
    IFFT = params.IFFT

    ys = pre.ys
    Fys = pre.Fys

    A = 1 / sqrt(params.dx)

    # Generate random phases
    randn!(rng, ys)

    @inbounds for i in eachindex(Fys)
        Fys[i] = ys[i]
    end

    FFT * Fys # In place FFT
    # n = length(Fys)
    # randn!(rng, @view Fys[1:n÷2+1])

    @inbounds for i in eachindex(Fys)
        Fys[i] *= A * d * sqrt(correlation(xks[i], params.surf))
    end

    # @inbounds for i in 1:n÷2+1
    #     Fys[i] *= A * d * sqrt(correlation(xks[i], params.surf))
    # end
    # # Enforce Hermitian symmetry for negative frequencies
    # @inbounds for i in (n÷2+2):n
    #     Fys[i] = conj(Fys[n-i+2])  # Conjugate symmetry F(-k) = F*(k)
    # end
    # # Special case for Nyquist frequency (if n is even)
    # if n % 2 == 0
    #     Fys[n÷2+1] = real(Fys[n÷2+1])
    # end
    IFFT * Fys # in place Inverse FFT

    @inbounds for i in eachindex(ys)
        ys[i] = real(Fys[i])
    end

    return nothing
end
function generate_surface!(pre::Preallocated, ::Parameters{FlatSurface,_MA,_MB})::Nothing where {_MA,_MB}
    pre.ys .= 0.0
    return nothing
end
function generate_surface!(pre::Preallocated, params::Parameters{SingleBumpSurface,_MA,_MB})::Nothing where {_MA,_MB}
    surf = params.surf
    δ = surf.d
    a = surf.a

    pre.ys .= δ * exp.(-0.5 * (params.xs / a) .^ 2)
    return nothing
end