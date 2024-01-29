using CairoMakie
using LaTeXStrings
CairoMakie.activate!(type="svg")
function read_ComplexF64_binary(fn, Nq)
    open(fn, "r") do io
        data = reshape(reinterpret(ComplexF64, read(io)), Nq, :)
        return data
    end
end

function read_Float64_binary(fn)
    open(fn, "r") do io
        data = reshape(reinterpret(Float64, read(io)), :)
        return data
    end
end
# Read from data folder and config.txt
#=
data_dir = readdir("data")
popat!(data_dir, findfirst(isequal("data-processed"), data_dir))
display(data_dir)

mdrcs_incoh = Vector{Vector{Float64}}(undef, length(data_dir))

for i in eachindex(data_dir)
    mdrcs_incoh[i] = read_Float64_binary(joinpath("data/", data_dir[i], "incoh.bin"))
end
=#

θs = range(-90, 90, length=length(cvec))
θ0 = 10.0

fig = Figure(fontsize=24)
ax = Axis(fig[1, 1], xlabel=L"$\theta_s$ [deg]", ylabel=L"Incoh. MDRC $\left[10^{-8}\right]$")
pad = 1.1
scale = 10^8
ylim = pad * maximum(cvec) * scale
lines!(ax,
    θs, scale .* cvec,
    color=:red,
    linestyle=:solid,
    linewidth=1)
# Incidence angle
lines!(ax,
    [θ0, θ0],
    [0, ylim],
    label=L"$\pm$%$θ0",
    color=:black,
    linestyle=:dot,
    linewidth=2)
lines!(ax,
    [-θ0, -θ0],
    [0, ylim],
    color=:black,
    linestyle=:dot,
    linewidth=2)

axislegend()
# tightlimits!(ax))
limits!(ax, -90, 90, 0, ylim)
display(fig)

save("plots/mdrc_incoh_ε-7.5+0.24i_p-type.pdf", fig)
