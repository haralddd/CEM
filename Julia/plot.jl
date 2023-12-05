using DelimitedFiles
using Plots
using Statistics


function plot_surf_tests()
    ims = range(1e-1, 1e-4, 50)

    # res_flat = readdlm("res_flat.csv", ',')
    # res_bump = readdlm("res_bump.csv", ',')
    # res_gaussian = readdlm("res_gaussian.csv", ',')

    plot(ims, res_flat, label="Flat", xlabel="Imaginary part of ε", ylabel="Unitarity")
    plot!(ims, res_bump, label="Single bump")
    plot!(ims, res_gaussian, label="Gaussian bump") |> display
    nothing
end

function read_complex_binary(fn, Nq)
    open(fn, "r") do io
        data = reshape(reinterpret(ComplexF64, read(io)), Nq, :)
        return data
    end
end

α(q, εμ) = sqrt(εμ - q^2)
α0(q) = sqrt(1.0 - q^2)

function MDRC_coh(A)
    res = Vector{Float64}(undef, size(A, 1))
    for i in axes(A, 2)
        res[i] = maximum(abs2.(A[:, i]))
    end
end

function MDRC_incoh(A, k, qis, qs, L)


    res = Vector{Float64}(undef, length(qis))
    for (i, qi) in enumerate(qis)
        res[i] = 1 / (2π * L) * α0(qs[qi]) / α0(k) * (mean(abs2.(A[qi, :])) - abs2(mean(A[qi, :])))
    end
    return res
end

function plot_gglass()
    Nq = 1025
    fn = "data/2023-11-21T18:44:01.786/gglass_θ0.0_νs_Nq1024_Nens1000.bin"
    A0s = read_complex_binary(fn, Nq)
    fn = "data/2023-11-21T18:44:01.786/gglass_θ34.05_νs_Nq1024_Nens1000.bin"
    A1s = read_complex_binary(fn, Nq)

    fn = "data/2023-11-23T11:36:20.385/gglass_θ0.0_νp_Nq1024_Nens1000.bin"
    A0p = read_complex_binary(fn, Nq)
    fn = "data/2023-11-23T11:36:20.385/gglass_θ34.05_νp_Nq1024_Nens1000.bin"
    A1p = read_complex_binary(fn, Nq)

    Q = 4
    Δq = Q / Nq
    qs = -Q/2:Δq:Q/2
    qis = findall(x -> -1.0 < x < 1.0, qs)
    c = 299792458.0 # Speed of light
    λ = 632.8e-9 # He-Ne laser wavelength
    L = 50 * λ # Length of surface
    ω = 2π * c / λ
    # L *= ω / c

    k0 = sind(0)
    k1 = sind(34.05)

    plt = plot(qs[qis], MDRC_incoh(A0s, k0, qis, qs, L), label="s-type 0°")
    plot!(qs[qis], MDRC_incoh(A0p, k0, qis, qs, L), label="p-type 0°")

    savefig(plt, "MDRC_incoh_s.pdf")

    plt = plot(qs[qis], MDRC_incoh(A1s, k1, qis, qs, L), label="s-type 34.05°")
    plot!(qs[qis], MDRC_incoh(A1p, k1, qis, qs, L), label="p-type 34.05°")

    savefig(plt, "MDRC_incoh_p.pdf")
    nothing
end

function plot_gsilver()
    Nq = 1025
    fn = "data/2023-11-23T15:02:13.270/gsilver_θ0.0_νs_Nq1024_Nens1000.bin"
    A0s = read_complex_binary(fn, Nq)
    fn = "data/2023-11-23T15:02:13.270/gsilver_θ34.05_νs_Nq1024_Nens1000.bin"
    A1s = read_complex_binary(fn, Nq)

    fn = "data/2023-11-23T14:59:05.366/gsilver_θ0.0_νp_Nq1024_Nens1000.bin"
    A0p = read_complex_binary(fn, Nq)
    fn = "data/2023-11-23T14:59:05.366/gsilver_θ34.05_νp_Nq1024_Nens1000.bin"
    A1p = read_complex_binary(fn, Nq)

    Q = 4
    Δq = Q / Nq
    qs = -Q/2:Δq:Q/2
    qis = findall(x -> -1.0 < x < 1.0, qs)
    c = 299792458.0 # Speed of light
    λ = 632.8e-9 # He-Ne laser wavelength
    L = 50 * λ # Length of surface
    ω = 2π * c / λ
    # L *= ω / c

    k0 = sind(0)
    k1 = sind(34.05)

    plt = plot(qs[qis], MDRC_incoh(A0s, k0, qis, qs, L), label="s-type 0°")
    plot!(qs[qis], MDRC_incoh(A0p, k0, qis, qs, L), label="p-type 0°")

    savefig(plt, "MDRC_incoh_s_silver.pdf")

    plt = plot(qs[qis], MDRC_incoh(A1s, k1, qis, qs, L), label="s-type 34.05°")
    plot!(qs[qis], MDRC_incoh(A1p, k1, qis, qs, L), label="p-type 34.05°")

    savefig(plt, "MDRC_incoh_p_silver.pdf")
    nothing
end

# plot_gglass()
plot_gsilver()