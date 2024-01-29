# using MKL
push!(LOAD_PATH, "Julia/RayleighSolver/")
using RayleighSolver
using Dates
using CairoMakie
using LaTeXStrings

function setup_dir(str="data")
    display("Adding suffix to dir")
    str *= "_" * Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")
    mkdir(str)
    display("New dir: \"$str\"")
    return str
end

function plot_mdrc(coh, incoh, θ0s=[0, 10, 20], name="simonsen_fig18_p")
    θs = range(-90, 90, length=size(incoh, 1))

    fig = Figure(fontsize=24)
    ax = Axis(fig[1, 1], xlabel=L"$\theta_s$ [deg]", ylabel=L"Incoh. MDRC $$")
    ylim = maximum(incoh)
    for i in axes(incoh, 2)
        θ0 = θ0s[i]
        y = incoh[:, i]
        lines!(ax,
            θs, y,
            label=L"$\theta_0=$%$θ0",
            linestyle=:solid,
            linewidth=1)
        lines!(ax,
            [θ0, θ0],
            [0, ylim],
            color=:black,
            linestyle=:dot,
            linewidth=2)
        lines!(ax,
            [-θ0, -θ0],
            [0, ylim],
            color=:black,
            linestyle=:dot,
            linewidth=2)

    end
    axislegend()
    limits!(ax, -90, 90, 0, ylim)
    display(fig)

    save("plots/$name.pdf", fig)
end

function cmd_line_run()
    println("Running as script...")
    print("Make new run configuration [y] or load existing [any]? ")
    in = readline()
    if in == "y"
        rp, sp, generator! = make_solver_config()
    else
        print("Enter config file name [./input/{config.txt}]: ")
        name = readline()
        rp, sp, generator! = load_solver_config("input/" * (name == "" ? "example.txt" : name))
    end

    print("Enter save data directory [./data/{simulation_name}] (=$name): ")
    in = readline()
    run_dir = "data/" * (in == "" ? name : in)
    run_dir = setup_dir(run_dir)
    save_solver_config(run_dir * "/config.txt", rp)

    print("Enter number of ensemble averages [100]: ")
    in = readline()
    N_ens = parse(Int, in == "" ? "100" : in)
    open(run_dir * "/N_ens.bin.int64", "w") do io
        write(io, N_ens)
    end
    coh, incoh = solve_MDRC!(rp, sp, generator!, N_ens)

    # Save data
    for i in axes(incoh, 2)
        open(run_dir * "/coh$i.bin", "w") do io
            write(io, coh)
        end
        open(run_dir * "/incoh$i.bin", "w") do io
            write(io, incoh)
        end
    end
    display("Saved data to $run_dir")
end


function predefined_run()
    display("Running for p")
    rp, sp, generator! = load_solver_config("input/simonsen_fig18_p.txt")

    coh, incoh = solve_MDRC!(rp, sp, generator!, 10000)
    display("Plotting for p")
    plot_mdrc(coh, incoh, [0, 10, 20], "simonsen_fig18_p")

    display("Running for s")
    rp, sp, generator! = load_solver_config("input/simonsen_fig18_s.txt")

    coh, incoh = solve_MDRC!(rp, sp, generator!, 10000)
    plot_mdrc(coh, incoh, [0, 10, 20], "simonsen_fig18_s")
    display("Done and saved plots")
end

if abspath(PROGRAM_FILE) == @__FILE__
    print("Load a config [y] or run as defined in code [any]? ")
    in = readline()
    if in == "y"
        cmd_line_run()
    else
        predefined_run()
    end
else
    predefined_run()
end

