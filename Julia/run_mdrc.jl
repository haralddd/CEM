# using MKL
push!(LOAD_PATH, "$(@__DIR__)/RayleighSolver/")
using RayleighSolver
using Dates
using CairoMakie
using LaTeXStrings
using Statistics

function timestamp_suffix(str="data")
    return str * "_" * Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")
end

function plot_mdrc(; coh, incoh, θ0s=[0, 10, 20], colors=[:red, :blue, :green], name="simonsen_fig18_p")
    θs = range(-90, 90, length=size(incoh, 1))

    fig = Figure(fontsize=24)
    ax = Axis(fig[1, 1], xlabel=L"$\theta_s$ [deg]", ylabel=L"Incoh. MDRC $$")
    # ymean = mean(incoh)
    # ymax = maximum(incoh)
    # ylim = 1.1ymax > 2.0 * ymean ? 1.1ymax : 2.0 * ymean
    for i in axes(incoh, 2)
        θ0 = θ0s[i]
        y = incoh[:, i]
        lines!(ax,
            θs, y,
            label=L"$\theta_0=$%$θ0",
            linestyle=:solid,
            linewidth=1,
            color=colors[i])
        vlines!(ax,
            [-θ0, θ0],
            color=colors[i],
            linestyle=:dot,
            linewidth=1)

    end
    axislegend()

    mkpath("plots")
    save("plots/$name.pdf", fig)
    println("Saved plot to plots/$name.pdf")
    display(fig)

end

function cmd_line_run()
    println("Running as script...")
    print("Make new run configuration [=y] or load existing [n]? ")
    input = readline()
    input = input == "" ? "y" : input
    if input == "y"
        spa = config_creation_prompt()
    else
        files = readdir("input")
        println("Choose config file to load [Int64] (=1): ")
        for i in eachindex(files)
            println("[$i] $(files[i])")
        end
        input = readline()
        input = input == "" ? "1" : input
        spa = load_spa_config("input/" * files[parse(Int64, input)])
    end
    print("Enter output label [string] (=\"default\"): ")
    input = readline()
    outlabel = input == "" ? "default" : input

    print("Do you want timestamp suffix? [y|n] (=y): ")
    input = readline()
    input = input == "" ? "y" : input
    outlabel = input == "y" ? timestamp_suffix(outlabel) : outlabel
    outpath = "data/" * outlabel * ".jld2"
    save_spa_config(outpath, spa)
    println("Saved config to $outpath")

    print("Enter number of ensembles to run (=100): ")
    input = readline()
    N_ens = parse(Int, input == "" ? "100" : input)
    save(outpath, "N_ens", N_ens)

    sp = SimPrealloc(spa)
    coh, incoh = solve_MDRC!(sp, spa, N_ens)
    save_mdrc_data(outpath, coh, incoh)
    println("Done and saved data to $outpath")

    print("Save plots? [y|n] (=y): ")
    input = readline()
    input = input == "" ? "y" : input
    if input == "y"
        plot_mdrc(; coh=coh, incoh=incoh, name=outlabel)
    end
end

function run_mdrc_silver()
    display("Running for p")
    spa, sp, generator! = load_solver_config("input/silver_rect_p")

    N = 10_000

    coh, incoh = solve_MDRC!(spa, sp, generator!, N)
    qis = findall(q -> q > -1 && q < 1, spa.qs)

    unit = zeros(size(coh, 2))
    for (i, qi) in enumerate(qis)
        unit[:] += (incoh[i, :] + coh[i, :]) / α0(spa.qs[qi])
    end
    display(unit)

    display("Plotting for p")
    plot_mdrc(; coh=coh, incoh=incoh, name="silver_rect_p")


    display("Running for s")
    spa, sp, generator! = load_solver_config("input/silver_rect_s")

    coh, incoh = solve_MDRC!(spa, sp, generator!, N)

    unit = zeros(size(coh, 2))
    for (i, qi) in enumerate(qis)
        unit[:] += (incoh[i, :] + coh[i, :]) / α0(spa.qs[qi])
    end
    display(unit)


    plot_mdrc(; coh=coh, incoh=incoh, name="silver_rect_s")
    display("Done and saved plots")
end

function run_mdrc_magnetic()
    N = 10_000

    display("Running for magnetic p gaussian")
    spa, sp, generator! = load_solver_config("input/magnetic_p")
    coh, incoh = solve_MDRC!(spa, sp, generator!, N)
    plot_mdrc(; coh=coh, incoh=incoh, name="magnetic_p_gaussian")


    display("Running for magnetic s gaussian")
    spa, sp, generator! = load_solver_config("input/magnetic_s")
    coh, incoh = solve_MDRC!(spa, sp, generator!, N)
    plot_mdrc(; coh=coh, incoh=incoh, name="magnetic_s_gaussian")

    display("Running for magnetic p rect")
    spa, sp, generator! = load_solver_config("input/magnetic_p_rect")
    coh, incoh = solve_MDRC!(spa, sp, generator!, N)
    plot_mdrc(; coh=coh, incoh=incoh, name="magnetic_p_rect")

    display("Running for magnetic s rect")
    spa, sp, generator! = load_solver_config("input/magnetic_s_rect")
    coh, incoh = solve_MDRC!(spa, sp, generator!, N)
    plot_mdrc(; coh=coh, incoh=incoh, name="magnetic_s_rect")

    display("Done and saved plots")
end

if abspath(PROGRAM_FILE) == @__FILE__
    print("Custom config [=y], run silver [silver], or run magnetic [magnetic]? ")

    input = readline()
    input = input == "" ? "y" : input
    if input == "y"
        cmd_line_run()
    elseif input == "silver"
        run_mdrc_silver()
    elseif input == "magnetic"
        run_mdrc_magnetic()
    else
        println("Invalid input")
    end
end

