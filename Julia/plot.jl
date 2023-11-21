using DelimitedFiles
using Plots


function plot_surf_tests()
    ims = range(1e-1, 1e-4, 50)

    res_flat = readdlm("res_flat.csv", ',')
    res_bump = readdlm("res_bump.csv", ',')
    res_gaussian = readdlm("res_gaussian.csv", ',')

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

function plot_gglass(dir)
    fn = "data/2023-11-21T18:44:01.786/gglass_θ0.0_νs_Nq1024_Nens1000.bin"
    A0s = read_complex_binary(fn, 1025)
    fn = "data/2023-11-21T18:44:01.786/gglass_θ0.0_νp_Nq1024_Nens1000.bin"
    A0p = read_complex_binary(fn, 1025)


end