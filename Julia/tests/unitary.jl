push!(LOAD_PATH, "Julia/RayleighSolver/")
using RayleighSolver
using CairoMakie
using LaTeXStrings

function unitary(R, rp, k)
    sum = 0.0
    idxs = findall(x -> -1.0 < x < 1.0, rp.qs)
    for i in idxs
        sum += abs2(R[i] * sqrt(α0(rp.qs[i]) / α0(k)))
    end
    return sum
end

function test_solvers(; surf_t::SurfType=flat, ε=2.25, μ=1.0, ν::Polarization=p)
    θ0 = 45.0
    L = 10.0e-6
    ims = range(1e-1im, 1e-4im, length=50)

    res = Vector{Float64}(undef, length(ims))
    for i in eachindex(ims)
        im = ims[i]
        rp = RayleighParams(
            ν=ν,
            Nq=2^10,
            ε=ε + (ν == p ? im : 0),
            μ=μ + (ν == s ? im : 0),
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


function plot_unitarity()
    ims_p, res_p = test_solvers(surf_t=flat, ε=-20, μ=1, ν=p)
    ims_s, res_s = test_solvers(surf_t=flat, ε=1, μ=-20, ν=s)

    fig = Figure(; size=(500, 400))
    ax = Axis(fig[1, 1], xlabel=L"Imaginary part of $\varepsilon$", ylabel=L"$$Unitarity")
    scatter!(ax, imag.(ims_p), res_p, label="p", color=(:blue, 0.5), marker=:circle)
    scatter!(ax, imag.(ims_s), res_s, label="s", color=(:red, 0.5), marker=:rect)
    axislegend(ax, loc=:tr)
    display(fig)
    save("plots/unitarity.pdf", fig)
end

plot_unitarity()
