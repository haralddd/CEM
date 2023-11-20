using MKL
push!(LOAD_PATH, "Julia/RayleighSetup/")
push!(LOAD_PATH, "Julia/RayleighSolver/")
using RayleighSetup
using RayleighSolver
using Plots
using BenchmarkTools

rp, sp_f, sp_g = test_rp_and_sp();

begin
    plot(rp.xs, sp_g.ys, linestyle=:solid, size=(400, 400), label="Gaussian surface")
    plot!(rp.xs, sp_f.ys, linestyle=:dash, label="Flat surface")
    plot!([rp.xs[begin], rp.xs[end]], [rp.δ, rp.δ], linestyle=:dot, linewidth=1, label="RMS height")

    xlabel!("x * ω/c [m]")

    ratio = 1e-1
    xlims!(-rp.Lx / 2, rp.Lx / 2)
    ylims!(-rp.Lx / 2 * ratio, rp.Lx / 2 * ratio) |> display
end




Nq = 2^11

silver = -17.5 + 0.48im
glass = 2.25

rp_p = RayleighParams(;
    ν=p,
    Nq=Nq,
    ε=glass,
    L=10.0e-6,
    Q_mult=4,
    Ni=10
);


reduction = x::Vector{ComplexF64} -> maximum(abs.(x))^2

# show_params(rp_p)
@time sol_p = solve(rp_p, reduction);

sol_p


# sum(abs2, sol_p)
rp_s = RayleighParams(;
    ν=s,
    Nq=Nq,
    ε=glass,
    L=10.0e-6,
    Q_mult=4,
    Ni=10
);
# show_params(rp_s)
@time sol_s = solve(rp_s);
θs_new_p = asind.(rp_p.ks)
θs_new_s = asind.(rp_s.ks)

refl_p = (maximum(abs.(sol_p), dims=1) |> vec) .^ 2
refl_s = (maximum(abs.(sol_s), dims=1) |> vec) .^ 2

display(refl_s[1])



ε = glass

r_analytical = abs2.((1.0 .- .√(ε .- rp_p.ks .^ 2)) ./ (1 .+ .√(ε .- rp_p.ks .^ 2)))

brewster = 56

bi = findfirst(refl_p .== minimum(refl_p))
display("Brewster angle $brewster vs. $(θs_new_p[bi])")

plot(θs_new_p, refl_p, label="ν = p", marker=(:circle, 2));
plot!(θs_new_s, refl_s, label="ν = s", marker=(:circle, 2));
plot!(θs_new_p, r_analytical, label="Analytical Fresnel")
plot!([0.0, 90.0], [1.0, 1.0], linestyle=:dash, linecolor=:black, linewidth=1, label=nothing);
plot!([90.0, 90.0], [0.0, 1.0], linestyle=:dash, linecolor=:black, linewidth=1, label=nothing);


xlims!(0, 90);
ylims!(0, 1.00);
xlabel!("Θ [°]");
ylabel!("Fresnel reflection coefficient")


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

function surf_distr()
    L = 10.0e-6
    width = 10.0e-9
    height = 50.0e-9
    rp = RayleighParams(
        ν=p,
        Nq=2^12,
        ε=-20.0 + 0.48im,
        Ni=10,
        L=L,
        δ=height,
        a=width,
    )
    N = 1000
    ys = Matrix{Float64}(undef, length(rp.xs), N)
    for n in 1:N
        ys[:, n] = SurfPreAlloc(rp, gaussian).ys
    end
    display(rp.δ / (rp.ω / c0))
    display("$(.√(mean((ys) .^ 2))) ≈ $(rp.δ) : $(.√(mean(ys .^ 2)) ≈ rp.δ)")
    # plt = plot(rp.xs, ys[:, 1], label=nothing, alpha=0.5)
    # for n in 2:N
    #     plot!(rp.xs, ys[:, n], label=nothing, alpha=0.5)
    # end
    # display(plt)
    # histogram(rp.xs, weights=mean(ys .^ 2; dims=2), bins=100, label="Gaussian surface", normalize=:false) |> display

end

surf_distr()
