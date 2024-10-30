# using MKL
push!(LOAD_PATH, "$(@__DIR__)/RayleighSolver/")
using RayleighSolver
using Dates
include("plot.jl")


function /(x::AbstractString, y::AbstractString)
    return joinpath(x, y)
end

const DEFAULT_PATH = "$(@__DIR__)/"
const DEFAULT_INPUT = DEFAULT_PATH / "input"
const DEFAULT_OUTPUT = DEFAULT_PATH / "output"
const DEFAULT_PLOTDIR = DEFAULT_PATH / "plots"

mkpath(DEFAULT_INPUT)
mkpath(DEFAULT_OUTPUT)
mkpath(DEFAULT_PLOTDIR)

function timestamp_suffix(str="out")
    return str * "_" * Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")
end

function args_findoption(option)::Union{Nothing,Int}
    return findfirst(isequal(option), ARGS)
end

function get_jld_files(path)
    return filter(x -> endswith(x, ".jld2"), readdir(path))
end

function interface_prompt()
    print("Interface type [flat|gaussian|singlebump|rect(West O'Donnell)] (=gaussian): ")
    input = readline()
    input = input == "" ? "gaussian" : input

    if input == "rect"
        print("West O'Donnell RMS height, δ [nm] (=5.0): ")
        input = readline()
        d = parse(Float64, input == "" ? "5.0" : input) * 1e-9

        print("West O'Donnell lower cutoff, k- [scaled to omega/c] (=0.8): ")
        input = readline()
        km = parse(Float64, input == "" ? "0.8" : input)

        print("West O'Donnell upper cutoff, k+ [scaled to omega/c] (=1.2): ")
        input = readline()
        kp = parse(Float64, input == "" ? "1.2" : input)

        return RectangularSurface(d, km, kp)
    elseif input == "gaussian"
        print("Gaussian RMS height, δ [nm] (=5.0): ")
        input = readline()
        d = parse(Float64, input == "" ? "5.0" : input) * 1e-9

        print("Gaussian correlation length, a [nm] (=1.0): ")
        input = readline()
        a = parse(Float64, input == "" ? "1.0" : input) * 1e-9

        return GaussianSurface(d, a)
    elseif input == "singlebump"
        print("Single bump RMS height, δ [nm] (=5.0): ")
        input = readline()
        d = parse(Float64, input == "" ? "5.0" : input) * 1e-9

        print("Single bump correlation length, a [nm] (=1.0): ")
        input = readline()
        a = parse(Float64, input == "" ? "1.0" : input) * 1e-9

        return SingleBumpSurface(d, a)
    else
        return FlatSurface()
    end
end

function material_prompt()
    print("Material type [vacuum|isotropic|uniaxialcrystal] (=vacuum): ")
    input = readline()
    input = input == "" ? "vacuum" : input
    if input == "uniaxialcrystal"
        print("Uniaxial crystal ε⟂ [complex] (=1.0): ")
        input = readline()
        eps_perp = parse(ComplexF64, input == "" ? "1.0" : input)

        print("Uniaxial crystal ε∥ [complex] (=1.0): ")
        input = readline()
        eps_para = parse(ComplexF64, input == "" ? "1.0" : input)

        print("Uniaxial crystal μ⟂ [complex] (=1.0): ")
        input = readline()
        mu_perp = parse(ComplexF64, input == "" ? "1.0" : input)

        print("Uniaxial crystal μ∥ [complex] (=1.0): ")
        input = readline()
        mu_para = parse(ComplexF64, input == "" ? "1.0" : input)

        return UniaxialCrystal(eps_perp, eps_para, mu_perp, mu_para)
    elseif input == "isotropic"
        print("Isotropic ε [complex] (=1.0): ")
        input = readline()
        eps = parse(ComplexF64, input == "" ? "1.0" : input)

        print("Isotropic μ [complex] (=1.0): ")
        input = readline()
        mu = parse(ComplexF64, input == "" ? "1.0" : input)

        return Isotropic(eps, mu)
    end

    return Vacuum()
end

function config_creation_prompt(path=DEFAULT_INPUT)::SimParams
    print("Input for solver parameters input SimParams struct\n")
    print("Polarization [p|s] (=p): ")
    input = readline()
    nu = parse(Polarization, input == "" ? "p" : input)

    print("lambda [nm] (=632.8): ")
    input = readline()
    lambda = parse(Float64, input == "" ? "632.8" : input) * 1e-9

    print("Q [multiple of omega/c] (=4): ")
    input = readline()
    Q = parse(Int64, input == "" ? "4" : input)

    print("Nq (=1024): ")
    input = readline()
    Nq = parse(Int64, input == "" ? "1024" : input)

    print("angles [list of deg \'0,1,2...\' OR range \'0:1:10\' OR \'fresnel\'] (=0:10:20): ")
    input = readline()
    if input == "fresnel"
        angles = 0.0:0.5:90.0
    elseif contains(input, ':') # Range input
        from, step, to = parse.(Float64, split(input == "" ? "0:10:20" : input, ':'))
        angles = from:step:to
    else
        angles = parse.(Float64, split(input == "" ? "0, 10, 20" : input, ','))
    end

    print("Lx [multiple of lambda] (=100): ")
    input = readline()
    Lx = parse(Float64, input == "" ? "100" : input) * lambda

    print("Ni (=10): ")
    input = readline()
    Ni = parse(Int64, input == "" ? "10" : input)

    surf = interface_prompt()

    println("Material above the interface (+)")
    above = material_prompt()
    println("Material below the interface (-)")
    below = material_prompt()

    print("Seed [Int64 >= 0] (= random seed): ")
    input = readline()
    seed = parse(Int64, input == "" ? "-1" : input)

    spa = SimParams{typeof(surf),typeof(nu),typeof(above),typeof(below)}(
        lambda=lambda,
        Q=Q,
        Nq=Nq,
        ks=sind.(angles),
        Lx=Lx,
        Ni=Ni,
        surf=surf,
        above=above,
        below=below,
        seed=seed,
        rescale=true,
    )

    print("Save config as input file? [y|n] (=n): ")
    input = readline()
    if input == "y"
        print("Filename [default.jld2]: ")
        input = readline()
        input = input == "" ? "default.jld2" : input
        save_spa_config(path / input, spa,
            override=Dict(:seed => seed) # Override the seed to generate random seed using the input
        )
    end

    return spa
end

function cli_help()
    println("""

    MDRC Commands:
    help - print this help
    new - create new run configuration
    input [path] - load configuration file or directory at \'path\', defaults to \'./input\'
    info [Int|name] - show information about configuration file
    run [Int|name] - run a configuration file
    list - show list of loadable configurations in \'path\'
        and optionally run or show information

    """)
    return
end

function cli_list(path)
    files = get_jld_files(path)
    for (idx, file) in enumerate(files)
        println("$idx: $file")
    end
    return
end

function cli_info(filepath)
    iters = 0
    try
        iters = load_ensemble_iters(filepath)
    catch
        print("No number of ensemble iterations found in file '$(filepath)', specify number:")
        iters = parse(Int64, readline())
        @assert iters > 0 "Number of ensemble iterations must be > 0"
        save_ensemble_iters(filepath, iters)
    end
    try
        spa = load_spa_config(filepath)
        display(spa)
        display("Ensemble iterations: $iters")
    catch
        error("File '$(filepath)' not found or file not valid")
    end
    return
end

function cli_run(filepath)
    cli_info(filepath)
    spa = load_spa_config(filepath)
    iters = load_ensemble_iters(filepath)

    @info "Initializing SolverData and precomputing..."
    data = SolverData(spa, iters)

    @info "Solving system of equations..."
    @time solve_MDRC!(data)

    print("Custom label (datetime is appended): ")
    input = readline()
    fname = timestamp_suffix(input)

    save_spa_config(DEFAULT_OUTPUT / fname, spa)
    save_ensemble_iters(DEFAULT_OUTPUT / fname, iters)
    save_mdrc_data(DEFAULT_OUTPUT / fname, data.out)
    
    @info "Config and results saved to '$(DEFAULT_OUTPUT / fname)'"
    include("plot.jl")
    plot_mdrc(data, fname, DEFAULT_PLOTDIR)

    @info "Plots saved to '$(DEFAULT_PLOTDIR / fname)'"
    return
end

function cli_file(arg, argval, confpath)
    file_name = ""
    try
        files = get_jld_files(confpath)
        file_name = files[parse(Int64, argval)]
    catch
        file_name = argval
    end

    filepath = confpath / file_name

    if (arg == "info")
        cli_info(filepath)
        exit(0)
    elseif (arg == "run")
        cli_run(filepath)
        exit(0)
    else
        error("Unknown single file argument: $(arg)")
        exit(0)
    end
end



function cli_directory(input_path)
    if (idx = args_findoption("list")) !== nothing
        cli_list(input_path)
        println("""
        Actions:
            info [Int] - show information about configuration file
            run [Int] - run configuration""")
        input = readline()

        arg, argval = split(input, " ")
        cli_file(arg, argval, input_path)
    else
        error("Unknown action following input path: $(input_path)")
        exit(0)
    end
end

function cli_main()
    @debug "ARGS: $ARGS" 
    if length(ARGS) < 1 || ARGS[1] == "help"
        cli_help()
        exit(0)
    end

    if (idx = args_findoption("new")) !== nothing
        config_creation_prompt()
        exit(0)
    end

    input_path = DEFAULT_INPUT

    if (idx = args_findoption("input")) !== nothing
        @assert length(ARGS) > idx "Missing argument after 'input'"
        try
            input_path = parse(String, ARGS[idx+1])
        catch
            error("Invalid positional value after 'input': $(ARGS[idx+1])")
            exit(0)
        end
    end

    @debug input_path
    @debug isdir(input_path)
    @debug pwd()

    if (idx = args_findoption("info") !== nothing)
        if isfile(input_path)
            cli_info(input_path)
        elseif isdir(input_path)
            @assert length(ARGS) > idx "Missing argument after 'info'"
            cli_file("info", ARGS[idx+1], input_path)
        end
    elseif (idx = args_findoption("run") !== nothing)
        if isfile(input_path)
            cli_run(input_path)
        elseif isdir(input_path)
            @assert length(ARGS) > idx "Missing argument after 'run'"
            cli_file("run", ARGS[idx+1], input_path)
        end
    elseif isdir(input_path)
        cli_directory(input_path)
    end
end