using CairoMakie

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
p_dir = "data/2024-01-22T13:42:05.401"
Nq = 1024 + 1 # Add 1 because of the zero frequency
s_dir = "2024-01-22T"

p_θ1_incoh_fp = joinpath(p_dir, "θ0.0_mdrc_incoh.bin")
p_θ1_incoh_mdrc = read_Float64_binary(p_θ1_incoh_fp)

θs = range(-90, 90, length=length(p_θ1_incoh_mdrc))
# λ = 632.8e-9 # He-Ne laser wavelength
# c0 = 299792458 # Speed of light in vacuum

fig = Figure(resolution=(600, 400))
ax = Axis(fig[1, 1], xlabel="θ [deg]", ylabel="MDRC [m⁻¹]")
lines!(ax,
    θs, p_θ1_incoh_mdrc,
    label="p-polarization",
    color=:blue)
display(fig)