

"O'Donnel rectangular correlation function for surface generation"
H(x::Float64)::Float64 = (x > 0.0 ? 1.0 : 0.0)
Wr(x::Float64, km::Float64, kp::Float64)::Float64 = (sin(kp * x) - sin(km * x)) / ((kp - km) * x) # Handle when x = 0 in rect_gen!
gr(k::Float64, km::Float64, kp::Float64)::Float64 = π / (kp - km) * (H(kp - k) * H(k - km) + H(kp - k) * H(-k - km))

"Gaussian correlation function for surface generation"
Wg(x::Float64, a::Float64)::Float64 = exp(-x^2 / a^2)
gg(k::Float64, a::Float64)::Float64 = √π * a * exp(-(a * k / 2.0)^2)


function generate_impl!(::Type{T}, sp::SimulationPreAlloc, rp::RayleighParams)::Nothing where {T<:SurfaceParams}
    @error "generate_impl! not implemented for $(typeof(T))"
    return nothing
end

function generate_impl!(::Type{FlatSurfaceParams}, sp::SimulationPreAlloc, rp::RayleighParams)::Nothing
    sp.ys .= 0.0
    return nothing
end

function generate_impl!(::Type{GaussianSurfaceParams}, sp::SimulationPreAlloc, rp::RayleighParams)::Nothing
    surf = rp.surf
    δ = surf.d
    a = surf.a
    sp.Z .= δ * randn(rp.rng, Float64, length(sp.ys)) |> complex
    sp.ys .= rp.FFT \ ((rp.FFT * sp.Z) .* sqrt.(rp.FFT * Wg.(rp.xs, a))) .|> real
    return nothing
end

function generate_impl!(::Type{SingleBumpSurfaceParams}, sp::SimulationPreAlloc, rp::RayleighParams)::Nothing
    surf = rp.surf
    δ = surf.d
    a = surf.a

    sp.ys .= δ * exp.((-0.5 / a^2) * rp.xs .^ 2)
    return nothing
end

function generate_impl!(::Type{RectSurfaceParams}, sp::SimulationPreAlloc, rp::RayleighParams)::Nothing
    surf = rp.surf
    δ = surf.d
    km = surf.km
    kp = surf.kp

    sp.Z .= δ * randn(rp.rng, Float64, length(sp.ys)) |> complex
    zr = ceil(Int64, length(rp.xs) / 2)

    corr = Vector{ComplexF64}(undef, length(rp.xs))
    corr[zr] = 1.0 # From Taylor series around x = 0
    corr[zr+1:end] .= gr.(rp.xs[zr+1:end], km, kp)
    corr[1:zr-1] .= gr.(rp.xs[1:zr-1], km, kp)

    sp.ys .= rp.FFT \ ((rp.FFT * sp.Z) .* sqrt.(rp.FFT * corr)) .|> real
    return nothing
end

# Surface generation dispatch
function generate!(sp::SimulationPreAlloc, rp::RayleighParams)::Nothing
    generate_impl!(typeof(rp.surf), sp, rp)
    return nothing
end