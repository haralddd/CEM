using MKL
push!(LOAD_PATH, "Julia/RayleighSetup/")
using RayleighSetup
using FFTW
using Plots
using BenchmarkTools
using Caching

using LinearAlgebra
BLAS.get_config() |> display

function α(q::Float64, εμ::ComplexF64)::ComplexF64
    # Calculate the α parameter for a given p
    return sqrt(εμ - q^2)
end
function α0(q::Float64)::ComplexF64
    # Calculate the α0 parameter for a given q
    return sqrt(complex(1.0 - q^2))
end

function calc_M_invariant!(M::Array{ComplexF64,3}, rp::RayleighParams)::Nothing
    # Calculate the surface invariant part of the Mpq matrix
    # Invariant under surface change), but NOT incident angle θ0

    ps = rp.ps
    qs = rp.qs
    κ = rp.ν == p ? rp.ε : rp.μ
    εμ = rp.ε * rp.μ

    @inbounds for I in eachindex(IndexCartesian(), M)
        i, j, n = Tuple(I)
        p = ps[i]
        q = qs[j]
        M[i, j, n] = (
            (-1.0im)^(n - 1) / factorial(n - 1) * (
                (p + κ * q) * (p - q) * (α(p, εμ) - α0(q))^(n - 2) +
                (α(p, εμ) + κ * α0(q)) * (α(p, εμ) - α0(q))^(n - 1))
        )
    end
    return nothing
end

function calc_N_invariant!(N::Matrix{ComplexF64}, rp::RayleighParams, k::Float64)::Nothing
    # Calculate the surface invariant part of the Npk vector
    # Invariant under surface change, but NOT incident angle θ0
    # Input k must be the nearest value in qs for correct results
    ps = rp.ps
    κ = rp.ν == p ? rp.ε : rp.μ
    εμ = rp.ε * rp.μ


    @inbounds for I in eachindex(IndexCartesian(), N)
        (i, n) = Tuple(I)
        p = ps[i]
        N[i, n] = (
            (-1.0im)^(n - 1) / factorial(n - 1) * (
                (p + κ * k) * (p - k) * (α(p, εμ) + α0(k))^(n - 2) +
                (α(p, εμ) - κ * α0(k)) * (α(p, εμ) + α0(k))^(n - 1))
        )
    end
    return nothing
end

function solve_step!(sp::SurfPreAlloc, rp::RayleighParams,
    M_pre::Array{ComplexF64,3}, N_pre::Matrix{ComplexF64},
    ki::Int)::Vector{ComplexF64}
    sp.Mpq .= 0.0
    sp.Npk .= 0.0

    display("Calculating Mpq and Npk")
    @inbounds for n in axes(M_pre, 3)
        sp.Fys .= sp.ys .^ (n - 1)
        rp.FT * sp.Fys # In place FFT
        fftshift!(sp.sFys, sp.Fys) # No way to do shift in place apparently

        @inbounds for I in eachindex(IndexCartesian(), sp.Mpq)
            i, j = Tuple(I)
            sp.Mpq[i, j] += M_pre[i, j, n] * sp.sFys[i+j-1]
        end

        @inbounds for i in eachindex(sp.Npk)
            sp.Npk[i] -= N_pre[i, n] * sp.sFys[i+ki-1]
        end
    end

    display("Solving for R")
    return sp.Mpq \ sp.Npk
end

using IterTools
using LoopVectorization

Tuple(CartesianIndices(rand(3, 4)))

function solve!(sp::SurfPreAlloc, rp::RayleighParams)
    # Solve for the reflection coefficient
    # This function is not used in the final implementation
    # since it is more efficient to solve for R at each step
    sp.Mpq .= 0.0
    sp.Npk .= 0.0

    M(p, q, n) = (
        (-1.0im)^(n - 1) / factorial(n - 1) * (
            (p + κ * q) * (p - q) * (α(p, εμ) - α0(q))^(n - 2) +
            (α(p, εμ) + κ * α0(q)) * (α(p, εμ) - α0(q))^(n - 1))
    )

    N(p, k, n) = (
        (-1.0im)^(n - 1) / factorial(n - 1) * (
            (p + κ * k) * (p - k) * (α(p, εμ) + α0(k))^(n - 2) +
            (α(p, εμ) - κ * α0(k)) * (α(p, εμ) + α0(k))^(n - 1))
    )
    vmapreduce

    for n in 1:rp.Ni+1
        sp.Fys .= sp.ys .^ (n - 1)
        rp.FT * sp.Fys # In place FFT
        fftshift!(sp.sFys, sp.Fys) # No way to do shift in place apparently

        vmapntt!(I -> M(Tuple(I)..., n) * sp.sFys[I], sp.Mpq, eachindex(IndexCartesian(), sp.Mpq))
        vmapntt!(N, sp.Npk, product(rp.ps, n))
    end


    display("Solving for R")
    return sp.Mpq \ sp.Npk
end

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

function drc(R::ComplexF64, θs::Float64, θ0::Float64, rp::RayleighParams)
    # Calculate the differential reflection coefficient
    # ω/c factored out in rescaling, degrees in radians
    return cos(θs)^2 * abs2(R) * rp.ω / (2π * rp.Lx * c0 * cos(θ0))
end

function drc_q_space(R::ComplexF64, q::Float64, θ0::Float64, rp::RayleighParams)
    # Calculate the differential reflection coefficient
    # degrees in radians
    return abs2(R) * sqrt(1 - q^2) / (2π * rp.Lx * (c0 / rp.ω) * cos(θ0))
end

@enum SimType perf mem
function test_unitary(Nq=2^13 - 1, ε=-1.0, mem::SimType=perf)
    # This function checks the unitary condition of the scattering vector
    # for flat and single bump gaussian surfaces
    surf_flat::SurfType = flat::SurfType
    surf_bump::SurfType = singlebump::SurfType

    θ0_in = 45.0
    Nq = 2^13 - 1
    ε = -1.0

    L = 10.0e-6
    height = 30.0e-8 # Height of the bump
    width = L / 10.0 # Width of the bump

    rp = RayleighSetup.RayleighParams(
        ν=p,
        Nq=Nq,
        ε=ε,
        δ=height,
        a=width,
        L=L,
        Q_mult=4,
    )

    # Generate and show surfaces
    display("Generating surface and steps for flat surf")
    @time sp_flat = SurfPreAlloc(rp, surf_flat)

    display("Preallocating surface and steps for single bump surf")
    @time sp_bump = SurfPreAlloc(rp, surf_bump)

    display("Plotting surfaces for viz")
    plt = plot(rp.xs, sp_bump.ys, label="Bump")
    plot!(rp.xs, sp_flat.ys, label="Flat")
    ylims!(-10.0 * rp.δ, 10.0 * rp.δ)
    display(plt)

    # Calc the invariant parts
    display("Calculating invariant parts of Mpk")
    Mpq_pre = Array{ComplexF64,3}(undef, length(rp.ps), length(rp.qs), rp.Ni + 1)
    @time calc_M_invariant!(Mpq_pre, rp)


    ki = searchsortedfirst(rp.qs, sind(θ0_in), rev=true)
    k = rp.qs[ki]

    display("Calculating invariant parts of Npk")
    Npk_pre = Matrix{ComplexF64}(undef, length(rp.ps), rp.Ni + 1)
    @time calc_N_invariant!(Npk_pre, rp, k)

    display("Solving flat surface")
    sol_flat = solve_step!(sp_flat, rp, Mpq_pre, Npk_pre, ki)


    display("Solving single bump surface")
    sol_bump = solve_step!(sp_bump, rp, Mpq_pre, Npk_pre, ki)




    display("Unitary conditions")
    ridx = findall(x -> -1.0 < x < 1.0, rp.qs) # reduced domain of scattered light
    θ0 = asin(k)

    drc_pf = [drc_q_space(sol_flat[i], rp.qs[i], θ0, rp) for i in ridx]
    drc_pb = [drc_q_space(sol_bump[i], rp.qs[i], θ0, rp) for i in ridx]

    unit(R, q) = (rp.ω / c0)^2 * rp.Δq * abs2(R) * α0(q) / (α0(k) * 2π * rp.Lx)
    unit_pf = [unit(sol_flat[i], rp.qs[i]) for i in ridx]
    unit_pb = [unit(sol_bump[i], rp.qs[i]) for i in ridx]

    display("Flat p: $(sum(unit_pf))")
    display("Bump p: $(sum(unit_pb))")

    display("Flat p: $(sum(drc_pf)*rp.Δq)")
    display("Bump p: $(sum(drc_pb)*rp.Δq)")

    # show_params(rp)

    return nothing
end

test_unitary();