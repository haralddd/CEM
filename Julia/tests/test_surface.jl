push!(LOAD_PATH, "$(@__DIR__)/../RayleighSolver/")
using RayleighSolver
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
    valleys = xs[argminima(ys)]
    return mean(diff(peaks)) / 2

    acc = 0.0
    if peaks[1] < valleys[1] # First value is a peak
        for i in eachindex(valleys)
            acc += valleys[i] - peaks[i]
        end
        if length(peaks) > length(valleys) # Odd number of intervals
            acc += peaks[end] - valleys[end]
        end
    else # First value is a valley
        for i in eachindex(peaks)
            acc += peaks[i] - valleys[i]
        end
        if length(valleys) > length(peaks) # Odd number of intervals
            acc += peaks[end] - valleys[end]
        end
    end

    return acc / (length(peaks) + length(valleys))
end

function test_surf(surf::T) where T <:RandomSurface
    M = 10000
    spa, sp = default_params_for_surface_testing(surf)

    dx = spa.dx
    xs = spa.xs
    d = spa.surf.d

    means = zeros(M)
    rms = zeros(M)
    slopes = zeros(M)
    dists = zeros(M)
    for i in 1:M
        generate!(sp, spa)
        means[i] = mean(sp.ys)
        rms[i] = sqrt(mean(sp.ys .^ 2))
        slopes[i] = mean_slope(sp.ys, dx)
        dists[i] = mean_peak_valley_dist(xs, sp.ys)
    end



    H1 = "\033[92m"
    H2 = "\033[32m"
    reset_color = "\033[m"
    println_header(text) = println("$H2 - $text - $reset_color")

    println("$H1--- $surf ---$reset_color")

    println()
    println_header("Ensemble mean, ⟨ζ(x)⟩:")
    println("Analytical = 0.0")
    println("Numerical = $(mean(means))")

    println()
    println_header("Ensemble mean slope, s:")
    println("Analytical = $(analytical_slope(spa.surf))")
    println("Numerical = $(mean(slopes))")

    println()
    println_header("Ensemble mean peak-valley distance, ⟨D⟩:")
    println("Analytical = $(analytical_dist(spa.surf))")
    println("Numerical = $(mean(dists))")

    println()
    println_header("Ensemble mean RMS height, δ:")
    println("δ = $(d)")
    println("Numerical  = $(mean(rms))")
end

test_gaussian() = test_surf(GaussianSurface(30.0e-9, 100.0e-9))
test_rect() = test_surf(RectangularSurface(30.0e-9, 0.82, 1.97))

