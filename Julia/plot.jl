using CairoMakie
using LaTeXStrings

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
data_dir = readdir("data")
popat!(data_dir, findfirst(isequal("data-processed"), data_dir))
display(data_dir)

mdrcs_incoh1 = Vector{Vector{Float64}}(undef, length(data_dir))
mdrcs_incoh2 = Vector{Vector{Float64}}(undef, length(data_dir))

for i in eachindex(data_dir)
    mdrcs_incoh1[i] = read_Float64_binary(joinpath("data/", data_dir[i], "θ0.0_mdrc_incoh.bin"))
    mdrcs_incoh2[i] = read_Float64_binary(joinpath("data/", data_dir[i], "θ34.05_mdrc_incoh.bin"))
end

function save_plot(mdrc_p, mdrc_s, θs, θ0, title)
    fig = Figure(; size=(500, 400))
    ax = Axis(fig[1, 1], xlabel=L"$\theta_s$ [deg]", ylabel=L"$$MDRC")
    lines!(ax,
        θs, mdrc_p,
        label="p-polarization",
        color=:blue,
        linestyle=:solid,
        linewidth=1)
    lines!(ax,
        θs, mdrc_s,
        label="s-polarization",
        color=:red,
        linestyle=:solid,
        linewidth=1)
    # Incidence angle
    lines!(ax,
        [θ0, θ0],
        [0, maximum([mdrc_p; mdrc_s])],
        label=L"\pm\theta_0",
        color=(:black, 0.5),
        linestyle=:dot,
        linewidth=1)
    lines!(ax,
        [-θ0, -θ0],
        [0, maximum([mdrc_p; mdrc_s])],
        color=(:black, 0.5),
        linestyle=:dot,
        linewidth=1)

    axislegend()
    tightlimits!(ax)
    display(fig)
    save("plots/mdrc_$title.pdf", fig)
end

save_plot(mdrcs_incoh1[1], mdrcs_incoh1[2], range(-90, 90, length=length(mdrcs_incoh1[1])), 0, "ε2.25_μ1.0_θ0.0")
save_plot(mdrcs_incoh2[1], mdrcs_incoh2[2], range(-90, 90, length=length(mdrcs_incoh2[1])), 34.05, "ε2.25_μ1.0_θ34.05")

save_plot(mdrcs_incoh1[3], mdrcs_incoh1[4], range(-90, 90, length=length(mdrcs_incoh1[1])), 0, "ε-17.5+0.48im_μ1.0_θ0.0")
save_plot(mdrcs_incoh2[3], mdrcs_incoh2[4], range(-90, 90, length=length(mdrcs_incoh2[1])), 34.05, "ε-17.5+0.48im_μ1.0_θ34.05")

save_plot(mdrcs_incoh1[5], mdrcs_incoh1[6], range(-90, 90, length=length(mdrcs_incoh1[1])), 0, "ε1.0_μ-5.0+1.0im_θ0.0")
save_plot(mdrcs_incoh2[5], mdrcs_incoh2[6], range(-90, 90, length=length(mdrcs_incoh2[1])), 34.05, "ε1.0_μ-5.0+1.0im_θ34.05")

save_plot(mdrcs_incoh1[7], mdrcs_incoh1[8], range(-90, 90, length=length(mdrcs_incoh1[1])), 0, "ε1.0+1.0im_μ-5.0_θ0.0")
save_plot(mdrcs_incoh2[7], mdrcs_incoh2[8], range(-90, 90, length=length(mdrcs_incoh2[1])), 34.05, "ε1.0+1.0im_μ-5.0_θ34.05")