# using MKL
push!(LOAD_PATH, "Julia/RayleighSolver/")
using RayleighSolver
using Dates

function setup_dir(str="data")
    if !isdir(str)
        mkdir(str)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 5
        println("Usage: julia solve.jl [p|s] [ε] [μ] [N] [θ].
        p|s:\tPolarization
        ε:\tRelative permittivity of the material
        μ:\tRelative permeability of the material
        N:\tNumber of surface realizations to solve for
        θ:\tIncident angle in degrees")
        exit(1)
    end


    pol = polarization_from_string(ARGS[1])
    ε = parse(ComplexF64, ARGS[2])
    μ = parse(ComplexF64, ARGS[3])
    N = parse(Int, ARGS[4])
    θ = parse(Float64, ARGS[5])
else
    # Mostly for testing
    pol = s
    ε = -7.5 + 0.24im
    μ = 1.0
    N = 10000
    θ = 20.0
end

λ = 457.9e-9 # Wavelength of incident light
L = 100 * λ # Length of surface
δ = 5.0e-9 # RMS height of surface
a = 100e-9 # Correlation length of surface

Nq = 2^10
Ni = 10

rp = RayleighParams(
    ν=pol,
    λ=λ,
    Nq=Nq,
    ε=ε,
    μ=μ,
    Ni=Ni,
    L=L,
    δ=δ,
    a=a,
)


# Directory setup and dump config to txt
timestamp = now() |> string
run_dir = "data/" * timestamp
setup_dir(run_dir)
open(run_dir * "/config.txt", "w") do io
    write(io, params_as_string(rp))
    write(io, "\nN_ens = $N")
    write(io, "\nθ = $θ")
end

coh, incoh = run_threaded(rp; N_ens=N, θ=θ)

# Save data
open(run_dir * "/coh.bin", "w") do io
    write(io, coh)
end
open(run_dir * "/incoh.bin", "w") do io
    write(io, incoh)
end