# using MKL
push!(LOAD_PATH, "Julia/RayleighSolver/")
using RayleighSolver
using Dates
using CairoMakie
using LaTeXStrings
using Statistics

function setup_dir(str="data")
    display("Adding suffix to dir")
    str *= "_" * Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")
    mkdir(str)
    display("New dir: \"$str\"")
    return str
end

function save_data(coh, incoh, path)
    # Save data
    for i in axes(incoh, 2)
        open(path * "/coh$i.bin", "w") do io
            write(io, coh)
        end
        open(path * "/incoh$i.bin", "w") do io
            write(io, incoh)
        end
    end
    display("Saved data to $path")
end

function plot_mdrc(; coh, incoh, θ0s=[0, 10, 20], colors=[:red, :blue, :green], name="simonsen_fig18_p")
    θs = range(-90, 90, length=size(incoh, 1))

    fig = Figure(fontsize=24)
    ax = Axis(fig[1, 1], xlabel=L"$\theta_s$ [deg]", ylabel=L"Incoh. MDRC $$")
    ymean = mean(incoh)
    ymax = maximum(incoh)
    ylim = 1.1ymax > 2.0 * ymean ? 1.1ymax : 2.0 * ymean
    for i in axes(incoh, 2)
        θ0 = θ0s[i]
        y = incoh[:, i]
        lines!(ax,
            θs, y,
            label=L"$\theta_0=$%$θ0",
            linestyle=:solid,
            linewidth=1,
            color=colors[i])
        lines!(ax,
            [θ0, θ0],
            [0, ylim],
            color=:black,
            linestyle=:dot,
            linewidth=1)
        lines!(ax,
            [-θ0, -θ0],
            [0, ylim],
            color=:black,
            linestyle=:dot,
            linewidth=1)

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

    save_data(coh, incoh, run_dir)
end

function run_mdrc_silver()
    display("Running for p")
    rp, sp, generator! = load_solver_config("input/silver_rect_p")

    N = 10_000

    coh, incoh = solve_MDRC!(rp, sp, generator!, N)
    qis = findall(q -> q > -1 && q < 1, rp.qs)

    unit = zeros(size(coh, 2))
    for (i, qi) in enumerate(qis)
        unit[:] += (incoh[i, :] + coh[i, :]) / α0(rp.qs[qi])
    end
    display(unit)

    display("Plotting for p")
    plot_mdrc(; coh=coh, incoh=incoh, name="silver_rect_p")


    display("Running for s")
    rp, sp, generator! = load_solver_config("input/silver_rect_s")

    coh, incoh = solve_MDRC!(rp, sp, generator!, N)

    unit = zeros(size(coh, 2))
    for (i, qi) in enumerate(qis)
        unit[:] += (incoh[i, :] + coh[i, :]) / α0(rp.qs[qi])
    end
    display(unit)


    plot_mdrc(; coh=coh, incoh=incoh, name="silver_rect_s")
    display("Done and saved plots")
end

function run_mdrc_magnetic()
    N = 10_000

    display("Running for magnetic p gaussian")
    rp, sp, generator! = load_solver_config("input/magnetic_p")
    coh, incoh = solve_MDRC!(rp, sp, generator!, N)
    plot_mdrc(; coh=coh, incoh=incoh, name="magnetic_p_gaussian")


    display("Running for magnetic s gaussian")
    rp, sp, generator! = load_solver_config("input/magnetic_s")
    coh, incoh = solve_MDRC!(rp, sp, generator!, N)
    plot_mdrc(; coh=coh, incoh=incoh, name="magnetic_s_gaussian")

    display("Running for magnetic p rect")
    rp, sp, generator! = load_solver_config("input/magnetic_p_rect")
    coh, incoh = solve_MDRC!(rp, sp, generator!, N)
    plot_mdrc(; coh=coh, incoh=incoh, name="magnetic_p_rect")

    display("Running for magnetic s rect")
    rp, sp, generator! = load_solver_config("input/magnetic_s_rect")
    coh, incoh = solve_MDRC!(rp, sp, generator!, N)
    plot_mdrc(; coh=coh, incoh=incoh, name="magnetic_s_rect")

    display("Done and saved plots")
end

if abspath(PROGRAM_FILE) == @__FILE__
    print("Load a config [y], run silver [silver], or run magnetic [magnetic]? ")
    in = readline()
    if in == "y"
        cmd_line_run()
    elseif in == "silver"
        run_mdrc_silver()
    elseif in == "magnetic"
        run_mdrc_magnetic()
    else
        println("Invalid input")
    end
end

