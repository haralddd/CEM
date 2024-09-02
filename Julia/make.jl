#=
For command line use
=#
using Plots

push!(LOAD_PATH, "$(@__DIR__)/RayleighSolver/")
using RayleighSolver

if (abspath(PROGRAM_FILE) == @__FILE__)
    rp, sp = config_creation_prompt()
else
    rp, sp = config_default_creation()
    save_to(rp, "input/default")
    display(load_rp_struct("input/default.jld2"))
    # display(load_rp_desc("input/default.jld2"))
end

N_ens = 1000
coh, incoh = solve_MDRC!(rp, sp, N_ens)
display(coh)
Qs = range(start = -rp.Q, stop = rp.Q, length = size(coh,1))
fig1 = plot(Qs, log.(coh), legend = false, title = "N = $N_ens")
fig2 = plot(Qs, log.(incoh), legend = false, title = "N = $N_ens")
plt = plot(fig1, fig2, layout = (2, 1), size = (600, 800))
display(fig1)
display(fig2)