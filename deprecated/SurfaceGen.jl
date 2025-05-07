module SurfaceGen
export Surface, SurfaceEnsemble, rms_slope, ensemble_mean, ensemble_variance

include("DiscreteDiff.jl")
using .DiscreteDiff

using FFTW
using Statistics
using StaticArrays


struct Surface
    #=
    Makes a realization of a surface in the ensemble
    and saves the profile as well as derivatives:
    ζ[N], ζ'[N], ζ''[N]

    TODO: Does dual numbers work for this, might not be nescessary though
    See DualNumbers.jl
    =#

    ζ::Vector{Float64}
    ζ_dot::Vector{Float64}
    ζ_ddot::Vector{Float64}
end

struct SurfaceEnsemble
    surfs::Vector{Surface}
    xs::Vector{Float64}
    δ::Float64 # RMS height of surface profile function
    a::Float64 # Autocorrelation length
    L::Float64
    Δx::Float64
end

function SurfaceEnsemble(L::Float64, δ::Float64, a::Float64, N::Int, M::Int)
    W(x) = exp(-x^2 / a^2)
    g(k) = √π * a * exp(-(a * k / 2.0)^2)

    Δx = L / N
    xs = collect(-L/2.0:Δx:(L/2.0-Δx))

    P = plan_rfft(xs) # Pre-planned FFTs for performance
    # ks = rfftfreq(N, 2π / Δx)

    ζs = [Vector{Float64}(undef, N) for i in 1:M]
    W_p = W.(xs)

    for i in eachindex(ζs)
        ζs[i] .= P \ ((P * randn(Float64, N)) .* (P * W_p)) .* δ * 2π / sqrt(N)
    end

    # Take the finite difference and make surfaces
    surfs = [Surface(ζs[i], d1o2(ζs[i], Δx), d2o2(ζs[i], Δx)) for i ∈ 1:M]

    SurfaceEnsemble(surfs, xs, δ, a, L, Δx)
end

function rms_slope(surf::Surface)
    sqrt(mean((surf.ζ_dot) .^ 2))
end

end # Module SurfaceGen