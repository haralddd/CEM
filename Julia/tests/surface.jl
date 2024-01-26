push!(LOAD_PATH, "Julia/RayleighSolver/")
using RayleighSolver
using Statistics
using CairoMakie
using LaTeXStrings
using FFTW


function d1o2(x, Δx)
    # Discrete first derivative, order O(h²)
    # https://en.wikipedia.org/wiki/Finite_difference_coefficient
    dx = similar(x)
    N = length(dx)
    @assert N > 2

    # Forward on start point to order O(h²)
    dx[1] = -1.5 * x[1] + 2.0 * x[2] - 0.5 * x[3]
    # Backward on end points to order O(h²)
    dx[end] = 1.5 * x[end] - 2.0 * x[end-1] + 0.5 * x[end-2]

    # Central on the rest to order O(h²)
    for i in 2:N-1
        dx[i] = -0.5 * x[i-1] + 0.5 * x[i+1]
    end

    return dx ./ Δx
end

function col_map(f, A)
    # Apply a function to each column of a matrix
    out = similar(A)
    @inbounds for i in axes(A, 2)
        out[:, i] = f(A[:, i])
    end
    return out
end

function mean_slope(ys::Matrix{Float64}, Δx::Float64)
    # Mean slope of the surface, s
    return sqrt.(mean(col_map(y -> d1o2(y, Δx) .^ 2, ys)))
end

function mean_slope_int(g::Function, ks::Vector{Float64}, rp::RayleighParams)
    return rp.δ * sqrt(rp.Δq / 2π * sum(ks .^ 2 .* g.(ks, rp.a)))

end

function mean_dist(g::Function, ks::Vector{Float64}, a::Float64)
    # Mean distance of the surface, ⟨D⟩
    return π * sqrt(sum(ks .^ 2 .* g.(ks, a)) / (sum(ks .^ 4 .* g.(ks, a))))
end

function surf_distr()
    L = 20.0e-6
    δ = 100.0e-9
    a = 100.0e-9
    Nq = 2^10
    ε = -20.0 + 0.48im
    Ni = 10
    rp1 = RayleighParams(
        ν=p,
        Nq=Nq,
        ε=ε,
        Ni=Ni,
        L=L,
        δ=δ,
        a=a,
    )
    rp2 = RayleighParams(
        ν=s,
        Nq=2Nq,
        ε=ε,
        Ni=Ni,
        L=3L,
        δ=5δ,
        a=7a,
    )
    N = 10000
    ys1 = Matrix{Float64}(undef, length(rp1.xs), N)
    ys2 = Matrix{Float64}(undef, length(rp2.xs), N)
    for n in 1:N
        ys1[:, n] = SurfPreAlloc(rp1, gaussian).ys
        ys2[:, n] = SurfPreAlloc(rp2, gaussian).ys
    end

    Qs1 = -rp1.Q:rp1.Δq:rp1.Q |> collect
    Qs2 = -rp2.Q:rp2.Δq:rp2.Q |> collect

    display("Surface function and spectrum:")
    g1 = gg.(Qs1, rp1.a)
    W1 = rp1.FT * Wg.(rp1.xs, rp1.a) |> fftshift .|> abs

    plot(Qs1, W1, label="𝔽{W}(k)")
    plot!(Qs1, g1, label="g(k)") |> display

    display("--- Mean slope, s ---")
    display("Analytical = $(√2 * rp1.δ / rp1.a)")
    display("Analytical scaled = $(√2 * rp2.δ / rp2.a)")
    display("Numerical = $(mean_slope(ys1, rp1.Δx))")
    display("Numerical scaled: = $(mean_slope(ys2, rp2.Δx))")
    display("Integral = $(mean_slope_int(gg, Qs1, rp1))")
    display("Integral scaled = $(mean_slope_int(gg, Qs2, rp2))")

    display("--- Mean peak-valley distance, ⟨D⟩ ---")
    display("Analytical = $(π/√6 * rp1.a)")
    display("Analytical scaled = $(π/√6 * rp2.a)")
    display("Numerical = $(mean_dist(gg, Qs1, rp1.a))")
    display("Numerical scaled = $(mean_dist(gg, Qs2, rp2.a))")

    display("--- RMS height ---")
    display("Input = $(rp1.δ)")
    display("Input scaled = $(rp2.δ)")
    display("Result = $(.√(mean(ys1 .^ 2)))")
    display("Result scaled = $(.√(mean(ys2 .^ 2)))")


end

surf_distr()