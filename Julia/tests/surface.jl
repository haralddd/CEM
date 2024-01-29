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

function mean_slope(ys::Vector{Float64}, Δx::Float64)
    # Mean slope of the surface, s
    return d1o2(ys, Δx) |> x -> sqrt(mean(x .^ 2))
end

function mean_slope_all(ys::Matrix{Float64}, Δx::Float64)
    # Mean slope of the surface, s
    return sqrt.(mean(col_map(y -> d1o2(y, Δx) .^ 2, ys)))
end

function mean_slope_int(g::Function, ks::Vector{Float64}, Δq, δ, a)
    return δ * sqrt(Δq / 2π * sum(ks .^ 2 .* g.(ks, a)))

end

function mean_dist(g::Function, ks::Vector{Float64}, a::Float64)
    # Mean distance of the surface, ⟨D⟩
    return π * sqrt(sum(ks .^ 2 .* g.(ks, a)) / (sum(ks .^ 4 .* g.(ks, a))))
end

function test_gaussian()
    L = 20.0e-6
    δ = 100.0e-9
    a = 100.0e-9
    Nq = 2^10
    ε = -20.0 + 0.48im
    Ni = 10

    δ1 = δ
    δ2 = δ * 3
    a1 = a
    a2 = a * 7
    L1 = L
    L2 = L * 3
    Nq1 = Nq
    Nq2 = Nq * 2

    surf1 = RayleighSolver.Surface(gaussian, [δ1, a1])
    surf2 = RayleighSolver.Surface(gaussian, [δ2, a2])

    rp1 = RayleighParams(
        ν=p,
        Nq=Nq1,
        ε=ε,
        Ni=Ni,
        L=L1,
        surf=surf1,
    )
    rp2 = RayleighParams(
        ν=s,
        Nq=Nq2,
        ε=ε,
        Ni=Ni,
        L=L2,
        surf=surf2,
    )
    N = 10000
    generate1! = surfgen_func_selector!(rp1)
    generate2! = surfgen_func_selector!(rp2)

    sp1 = SurfPreAlloc(rp1)
    sp2 = SurfPreAlloc(rp2)

    rms1 = zeros(N)
    rms2 = zeros(N)

    slope1 = zeros(N)
    slope2 = zeros(N)
    for i in 1:N
        generate1!(sp1.ys)
        generate2!(sp2.ys)
        rms1[i] = sqrt(mean(sp1.ys .^ 2))
        rms2[i] = sqrt(mean(sp2.ys .^ 2))

        slope1[i] = mean_slope(sp1.ys, rp1.Δx)
        slope2[i] = mean_slope(sp2.ys, rp2.Δx)
    end

    Qs1 = -rp1.Q:rp1.Δq:rp1.Q |> collect
    Qs2 = -rp2.Q:rp2.Δq:rp2.Q |> collect

    #=
    display("Surface function and spectrum:")
    a1s = rp1.surf.surf_params[2]
    a2s = rp2.surf.surf_params[2]
    δ1s = rp1.surf.surf_params[1]
    δ2s = rp2.surf.surf_params[1]
    g1 = gg.(Qs1, a1s)
    W1 = rp1.FT * Wg.(rp1.xs, a1s) |> fftshift .|> abs

    plot(Qs1, W1, label="𝔽{W}(k)")
    plot!(Qs1, g1, label="g(k)") |> display
    =#


    display("--- Mean slope, s ---")
    display("Analytical = $(√2 * δ1 / a1)")
    display("Analytical scaled = $(√2 * δ2 / a2)")
    display("Ens Numerical = $(mean(slope1))")
    display("Ens Numerical scaled: = $(mean(slope2))")

    #=
    display("--- Mean peak-valley distance, ⟨D⟩ ---")
    display("Analytical = $(π/√6 * a1)")
    display("Analytical scaled = $(π/√6 * a2)")
    display("Numerical = $(mean_dist(gg, Qs1, a1))")
    display("Numerical scaled = $(mean_dist(gg, Qs2, a2))")
    =#

    scale = rp1.ω / c0
    # scale = 1
    display("--- RMS height ---")
    display("Input = $(δ1)")
    display("Input scaled = $(δ2)")
    display("Result = $(mean(rms1) / scale)")
    display("Result scaled = $(mean(rms2) / scale)")
end

test_gaussian()

function test_rect()
    L = 20.0e-6
    δ = 100.0e-9
    Nq = 2^10
    ε = -20.0 + 0.48im
    Ni = 10

    δ1 = δ
    δ2 = δ * 3
    km1 = 0.8
    km2 = 0.7
    kp1 = 1.2
    kp2 = 1.3
    L1 = L
    L2 = L * 3
    Nq1 = Nq
    Nq2 = Nq * 2

    surf1 = RayleighSolver.Surface(rect, [δ1, km1, kp1])
    surf2 = RayleighSolver.Surface(rect, [δ2, km2, kp2])

    rp1 = RayleighParams(
        ν=p,
        Nq=Nq1,
        ε=ε,
        Ni=Ni,
        L=L1,
        surf=surf1,
    )
    rp2 = RayleighParams(
        ν=s,
        Nq=Nq2,
        ε=ε,
        Ni=Ni,
        L=L2,
        surf=surf2,
    )
    N = 1
    generate1! = surfgen_func_selector!(rp1)
    generate2! = surfgen_func_selector!(rp2)

    sp1 = SurfPreAlloc(rp1)
    sp2 = SurfPreAlloc(rp2)

    rms1 = zeros(N)
    rms2 = zeros(N)

    slope1 = zeros(N)
    slope2 = zeros(N)

    for i in 1:N
        generate1!(sp1.ys)
        generate2!(sp2.ys)

        rms1[i] = sqrt(mean(sp1.ys .^ 2))
        rms2[i] = sqrt(mean(sp2.ys .^ 2))

        slope1[i] = mean_slope(sp1.ys, rp1.Δx)
        slope2[i] = mean_slope(sp2.ys, rp2.Δx)
    end



    display("--- Mean slope, s ---")

    slope = (δ, km, kp) -> δ / √3 * √(kp^2 + kp * km + km^2)

    scale = rp1.ω / c0
    display("Analytical = $(slope(δ1, km1, kp1))")
    display("Analytical scaled = $(slope(δ2, km2, kp2))")
    display("Ens Numerical = $(mean(slope1) / scale)")
    display("Ens Numerical scaled: = $(mean(slope2) / scale)")

    # scale = 1
    display("--- RMS height ---")
    display("Input = $(δ1)")
    display("Input scaled = $(δ2)")
    display("Result = $(mean(rms1) / scale)")
    display("Result scaled = $(mean(rms2) / scale)")
end


test_rect()