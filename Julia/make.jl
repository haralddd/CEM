#=
For command line use
=#
using Plots
using FFTW

push!(LOAD_PATH, "$(@__DIR__)/RayleighSolver/")
using RayleighSolver

if (abspath(PROGRAM_FILE) == @__FILE__)
    spa, sp = config_creation_prompt()
else
    spa, sp = config_default_creation()
    save_to(spa, "input/default")
    display(load_rp_struct("input/default.jld2"))
    # display(load_rp_desc("input/default.jld2"))
end