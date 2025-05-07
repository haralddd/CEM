include("testconfig.jl")
using LaTeXStrings
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


function rms(xs, mean_xs)
    return sqrt(mean((xs .- mean_xs) .^ 2))
end
function rms_slope(ys::Vector{Float64}, Δx::Float64)
    # RMS slope of the surface, s
    d1 = d1o2(ys, Δx)
    d1_mean = mean(d1)
    return rms(d1, d1_mean)
end

function mean_peak_valley_dist(xs::Vector{Float64}, ys::Vector{Float64})
    peaks = xs[argmaxima(ys)]
    return mean(diff(peaks)) / 2
end

function test_surf(data::SolverData)


    params = data.params
    alloc = Preallocated(data)
    M = data.iters
    surf = params.surf

    xs = params.xs
    dx = params.dx
    d = params.surf.d

    meanval = 0.0
    rmsval = 0.0
    slope = 0.0
    dist = 0.0
    for m in 1:M
        generate_surface!(alloc, params)
        cur_mean = mean(alloc.ys)
        meanval = observe(meanval, cur_mean, m)
        rmsval = observe(rmsval, rms(alloc.ys, cur_mean), m)
        slope = observe(slope, rms_slope(alloc.ys, dx), m)
        dist = observe(dist, mean_peak_valley_dist(xs, alloc.ys), m)
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
    println("Analytical = $(analytical_slope(params.surf))")
    println("Numerical = $(slope)")

    println()
    println_header("Ensemble mean peak-valley distance, ⟨D⟩:")
    println("Analytical = $(analytical_dist(params.surf))")
    println("Numerical = $(dist)")

    println()
    println_header("Ensemble mean RMS height, δ:")
    println("δ = $(d)")
    println("Numerical  = $(rmsval)")
end
function make_plot(data::SolverData, label="")
    params = data.params
    alloc = Preallocated(data)

    generate_surface!(alloc, params)

    yscale = 4

    k0 = 2π / params.lambda
    xs = params.xs ./ k0 .* 1e6
    ys = alloc.ys ./ k0 .* 1e6
    maxy = max(extrema(ys)...)
    ylims = (-yscale * maxy, yscale * maxy)
    plt = plot(xs, ys, xlabel=L"x_1\ [\mu m]", ylabel=L"\zeta(x_1)\ [\mu m]", ylims=ylims, label=nothing)
    savefig(plt, "plots/$label.pdf")
end

if (abspath(PROGRAM_FILE) == @__FILE__) 
    iters = 100000
    test_gaussian1() = test_surf(default_gaussian_config(iters, 200))
    test_gaussian2() = test_surf(default_gaussian_config(iters, 100))
    test_gaussian3() = test_surf(default_gaussian_config(iters, 50))
    test_rect() = test_surf(default_rectangular_config(iters))

    test_gaussian1()
    test_gaussian2()
    test_gaussian3()
    test_rect()

    using Plots
    using LaTeXStrings
    make_plot(SolverData(ParametersConfig(surf=GaussianSurface(30.0e-9, 400.0e-9))), "gaussian_surf")
    make_plot(SolverData(ParametersConfig(surf=RectangularSurface(30.0e-9, 0.10, 1.10))), "rect_surf")
end