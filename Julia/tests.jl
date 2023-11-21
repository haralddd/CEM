using MKL
push!(LOAD_PATH, "Julia/RayleighSetup/")
push!(LOAD_PATH, "Julia/RayleighSolver/")
using RayleighSetup
using RayleighSolver
using Statistics
using Plots
using BenchmarkTools

function test_fresnel(; ε=2.25)
    surf_t::SurfType = flat::SurfType

    θs = 0.0:1.0:90.0
    Nq = 2^10

    rp_p = RayleighSetup.RayleighParams(
        ν=p,
        Nq=Nq,
        ε=ε
    )

    rp_s = RayleighParams(
        ν=s,
        Nq=Nq,
        ε=ε
    )

    # Calc the invariant part of Mpk
    display("Calculating invariant parts of Mpk")
    Mpk_p_pre = calc_Mpq_invariant(rp_p)
    Mpk_s_pre = calc_Mpq_invariant(rp_s)

    kis = [searchsortedfirst(rp_p.qs, sind(θs[i]), rev=true) for i in eachindex(θs)] |> unique
    ks = rp_p.qs[kis]

    # Results in Fresnel coefficients
    rs_p = Vector{Float64}(undef, length(ks))
    rs_s = Vector{Float64}(undef, length(ks))

    sp_p = SurfPreAlloc(rp_p, surf_t)
    sp_s = SurfPreAlloc(rp_s, surf_t)

    for i in eachindex(ks)
        Npk_p_pre = calc_Npk_invariant(rp_p, ks[i])
        Npk_s_pre = calc_Npk_invariant(rp_s, ks[i])
        # Solve and insert the specular reflection coefficient
        rs_p[i] = solve_step!(sp_p, rp_p, Mpk_p_pre, Npk_p_pre, kis[i]) .|> abs |> maximum
        rs_s[i] = solve_step!(sp_s, rp_s, Mpk_s_pre, Npk_s_pre, kis[i]) .|> abs |> maximum
    end


    display("Plotting Fresnel coefficients")
    plot(θs, rs_p .^ 2, label="p", marker=(:circle, 2))
    plot!(θs, rs_s .^ 2, label="s", marker=(:square, 2))

    return nothing
end

function trapz(xs, ys)
    # Simple trapezoidal integration with variable step size
    @assert length(xs) == length(ys)
    res = 0.0

    for i in eachindex(xs)
        i == 1 && continue
        res += 0.5 * abs(xs[i] - xs[i-1]) * (ys[i] + ys[i-1])
    end
    return res
end

function unitary(R, rp, k)
    sum = 0.0
    idxs = findall(x -> -1.0 < x < 1.0, rp.qs)
    for i in idxs
        sum += abs2(R[i]) * real(α0(rp.qs[i])) / real(α0(k))
    end
    return sum
end

function test_solvers(surf_t::SurfType=flat)
    θ0 = 45.0
    L = 10.0e-6
    ims = range(1e-1im, 1e-4im, length=50)

    res = Vector{Float64}(undef, length(ims))
    Threads.@threads for i in eachindex(ims)
        im = ims[i]
        rp = RayleighParams(
            ν=p,
            Nq=2^10,
            ε=-20.0 + im,
            Ni=10,
            L=L,
        )
        sp = SurfPreAlloc(rp, surf_t)
        ki = searchsortedfirst(rp.qs, sind(θ0), rev=true)
        k = rp.qs[ki]

        # Pre-calculate the invariant parts of the M and N matrices
        M_pre = Array{ComplexF64,3}(undef, length(rp.ps), length(rp.qs), rp.Ni + 1)
        N_pre = Matrix{ComplexF64}(undef, length(rp.ps), rp.Ni + 1)
        pre_M_invariant!(M_pre, rp)
        pre_N_invariant!(N_pre, rp, k)

        solve_pre!(sp, rp, M_pre, N_pre, ki)
        res[i] = unitary(sp.R, rp, k)
    end
    return ims, res
end

using DelimitedFiles

function run_unitary_tests()
    test = rand(100)
    writedlm("test.csv", test, ',')

    ims, res_flat = test_solvers(flat)
    writedlm("res_flat.csv", res_flat, ',')
    _, res_bump = test_solvers(singlebump)
    writedlm("res_bump.csv", res_bump, ',')
    _, res_gaussian = test_solvers(gaussian)
    writedlm("res_gaussian.csv", res_bump, ',')
end

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

using FFTW
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
