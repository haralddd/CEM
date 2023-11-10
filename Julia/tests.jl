setup_path = "Julia/RayleighSetup/"

if setup_path ∉ Base.load_path()
    push!(LOAD_PATH, setup_path)
end

using RayleighSetup
using Plots

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