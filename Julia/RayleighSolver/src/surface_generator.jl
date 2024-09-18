

"O'Donnel rectangular correlation function for surface generation"
H(x::Float64)::Float64 = (x > 0.0 ? 1.0 : 0.0)
Wr(x::Float64, km::Float64, kp::Float64)::Float64 = (sin(kp * x) - sin(km * x)) / ((kp - km) * x) # Handle when x = 0 in rect_gen!

"Gaussian correlation function for surface generation"
Wg(x::Float64, a::Float64)::Float64 = exp(-(x / a)^2)
gg(k::Float64, a::Float64)::Float64 = √π * a * exp(-(a * k / 2.0)^2)

function correlation(k::Float64, surf::GaussianSurfaceParams)::Float64
    a = surf.a
    return √π * a * exp(-(a * k / 2.0)^2)
end
function correlation(k::Float64, surf::RectSurfaceParams)::Float64
    kp = surf.kp
    km = surf.km
    return π / (kp - km) * (H(kp - k) * H(k - km) + H(kp + k) * H(-k - km))
end

"""
Generic Fourier filtering function, requires correlation(k, rp.surf::T) to be implemented
"""
function generate_impl!(::Type{T}, sp::SimulationPreAlloc, rp::RayleighParams)::Nothing where T<:SurfaceParams
    d = rp.surf.d
    ks = rp.xks

    # Generate random phases
    randn!(rp.rng, sp.ys)
    sp.Z .= d * sp.ys
    rp.FFT * sp.Z # In place FFT

    # Square root of correlation
    sp.ys .= sqrt.(correlation.(ks, Ref(rp.surf)) / rp.dx)
    sp.Z .= sp.Z .* sp.ys # Filter
    rp.FFT \ sp.Z # Inverse FFT
    sp.ys .= real.(sp.Z)

    return nothing
end

function generate_impl!(::Type{FlatSurfaceParams}, sp::SimulationPreAlloc, rp::RayleighParams)::Nothing
    sp.ys .= 0.0
    return nothing
end

function generate_impl!(::Type{SingleBumpSurfaceParams}, sp::SimulationPreAlloc, rp::RayleighParams)::Nothing
    surf = rp.surf
    δ = surf.d
    a = surf.a

    sp.ys .= δ * exp.((-0.5 * rp.xs / a) .^ 2)
    return nothing
end


"""
Generates an instance of the surface of given type, defined in `rp::RayleighParams`
Write the generated surface directly to the sp.ys array

"""
function generate!(sp::SimulationPreAlloc, rp::RayleighParams)::Nothing
    generate_impl!(typeof(rp.surf), sp, rp)
    return nothing
end