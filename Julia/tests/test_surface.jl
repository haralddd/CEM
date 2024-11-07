include("testconfig.jl")
using Statistics
using Peaks

analytical_slope(::T) where T<:RandomSurface = @error "Analytical slope not implemented for $T"
analytical_slope(surf::GaussianSurface) = √2 * surf.d / surf.a
function analytical_slope(surf::RectangularSurface)
    d = surf.d
    kp = surf.kp
    km = surf.km
    return d / √3 * √(kp^2 + kp*km + km^2)
end

analytical_dist(::T) where T<:RandomSurface = @error "Analytical distance not implemented for $T"
analytical_dist(surf::GaussianSurface) = π / √6 * surf.a
function analytical_dist(surf::RectangularSurface)
    kp = surf.kp
    km = surf.km
    return π * √(5/3 * (kp^3 - km^3) / (kp^5 - km^5) )
end
function d1o2(y, Δx)
    # Discrete first derivative, order O(h²)
    # https://en.wikipedia.org/wiki/Finite_difference_coefficient
    dy = similar(y)
    N = length(dy)
    @assert N > 2

    # Forward on start point to order O(h²)
    dy[1] = -1.5 * y[1] + 2.0 * y[2] - 0.5 * y[3]
    # Backward on end points to order O(h²)
    dy[end] = 1.5 * y[end] - 2.0 * y[end-1] + 0.5 * y[end-2]

    # Central on the rest to order O(h²)
    for i in 2:N-1
        dy[i] = -0.5 * y[i-1] + 0.5 * y[i+1]
    end

    return dy ./ Δx
end

function col_map(f, A)
    # Apply a function to each column of a matrix
    out = similar(A)
    @inbounds for i in axes(A, 2)
        out[:, i] = f(A[:, i])
    end
    return out
end

function mean_slope(ys::Vector{Float64}, Δx::Float64)
    # Mean slope of the surface, s
    return d1o2(ys, Δx) |> x -> sqrt(mean(x .^ 2))
end

function mean_peak_valley_dist(xs::Vector{Float64}, ys::Vector{Float64})
    peaks = xs[argmaxima(ys)]
    return mean(diff(peaks)) / 2
end

function test_surf(data::SolverData)


    spa = data.spa
    sp = data.sp
    M = data.iters
    surf = spa.surf

    xs = spa.xs
    dx = spa.dx
    d = spa.surf.d

    meanval = 0.0
    rms = 0.0
    slope = 0.0
    dist = 0.0
    for m in 1:M
        generate_surface!(sp, spa)
        meanval = observe(meanval, mean(sp.ys), m)
        rms = observe(rms, sqrt(mean(sp.ys .^ 2)), m)
        slope = observe(slope, mean_slope(sp.ys, dx), m)
        dist = observe(dist, mean_peak_valley_dist(xs, sp.ys), m)
    end

    H1 = "\033[92m"
    H2 = "\033[32m"
    reset_color = "\033[m"
    println_header(text) = println("$H2 - $text - $reset_color")

    println("$H1--- $surf ---$reset_color")

    println()
    println_header("Ensemble mean, ⟨ζ(x)⟩:")
    println("Analytical = 0.0")
    println("Numerical = $(meanval)")

    println()
    println_header("Ensemble mean slope, s:")
    println("Analytical = $(analytical_slope(spa.surf))")
    println("Numerical = $(slope)")

    println()
    println_header("Ensemble mean peak-valley distance, ⟨D⟩:")
    println("Analytical = $(analytical_dist(spa.surf))")
    println("Numerical = $(dist)")

    println()
    println_header("Ensemble mean RMS height, δ:")
    println("δ = $(d)")
    println("Numerical  = $(rms)")
end
const iters = 100000
test_gaussian() = test_surf(default_gaussian_config(iters))
test_gaussian2() = test_surf(default_gaussian_config(iters, 200*632.8e-9))
test_gaussian3() = test_surf(default_gaussian_config(iters, 50*632.8e-9))
test_rect() = test_surf(default_rectangular_config(iters))

