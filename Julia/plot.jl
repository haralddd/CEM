using DelimitedFiles
using Plots

ims = range(1e-1, 1e-4, 50)

res_flat = readdlm("res_flat.csv", ',')
res_bump = readdlm("res_bump.csv", ',')
res_gaussian = readdlm("res_gaussian.csv", ',')

plot(ims, res_flat, label="Flat", xlabel="Imaginary part of ε", ylabel="Unitarity")
plot!(ims, res_bump, label="Single bump")
plot!(ims, res_gaussian, label="Gaussian bump") |> display