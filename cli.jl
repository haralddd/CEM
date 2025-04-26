# using MKL
const DEFAULT_PATH = "$(@__DIR__)"
import Base./
function /(x::AbstractString, y::AbstractString)
    return joinpath(x, y)
end

using RayleighSolver
using Dates

const DEFAULT_INPUT = DEFAULT_PATH / "input"
const DEFAULT_OUTPUT = DEFAULT_PATH / "output"

mkpath(DEFAULT_INPUT)
mkpath(DEFAULT_OUTPUT)

function timestamp_suffix(str="out")
    return str * "_" * Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")
end

function args_findoption(option)::Union{Nothing,Int}
    return findfirst(isequal(option), ARGS)
end

function get_files(path)
    return filter(x -> endswith(x, ".conf") || endswith(x, ".jld2"), readdir(path))
end

function interface_prompt()
    print("Interface type [1:flat|2:gaussian|3:singlebump|4:rect(West O'Donnell)] (=gaussian): ")
    input = readline()
    input = input == "" ? "2" : input

    if input == "4"
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
    elseif input == "2"
        print("Gaussian RMS height, δ [nm] (=5.0): ")
        input = readline()
        d = parse(Float64, input == "" ? "5.0" : input) * 1e-9

        print("Gaussian correlation length, a [nm] (=100.0): ")
        input = readline()
        a = parse(Float64, input == "" ? "100.0" : input) * 1e-9

        return GaussianSurface(d, a)
    elseif input == "3"
        print("Single bump RMS height, δ [nm] (=5.0): ")
        input = readline()
        d = parse(Float64, input == "" ? "5.0" : input) * 1e-9

        print("Single bump correlation length, a [nm] (=100.0): ")
        input = readline()
        a = parse(Float64, input == "" ? "100.0" : input) * 1e-9

        return SingleBumpSurface(d, a)
    else
        return FlatSurface()
    end
end

function material_prompt(;default = "vacuum")
    print("Material type [1:vacuum|2:isotropic|3:uniaxial] (default: $default): ")
    input = readline()
    input = input == "" ? default : input
    if input == "3" || input == "uniaxial"
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

        return Uniaxial(eps_perp, eps_para, mu_perp, mu_para)
    elseif input == "2" || input == "isotropic"
        print("Isotropic ε [complex] (=1.0): ")
        input = readline()
        eps = parse(ComplexF64, input == "" ? "1.0" : input)

        print("Isotropic μ [complex] (=1.0): ")
        input = readline()
        mu = parse(ComplexF64, input == "" ? "1.0" : input)

        return Isotropic(eps, mu)
    elseif input == "1" || input == "vacuum"
        return Vacuum()
    else
        error("Unknown material type: $(input)")
    end
end

function config_creation_prompt(path=DEFAULT_INPUT)::Nothing
    if args_findoption("default") !== nothing
        save_parameters_config(path / "default", ParametersConfig())
        return
    end
    print("Input for solver parameters input Parameters struct\n")

    print("lambda [nm] (=632.8): ")
    input = readline()
    lambda = parse(Float64, input == "" ? "632.8" : input) * 1e-9

    print("Nx (=2048): ")
    input = readline()
    Nx = parse(Int64, input == "" ? "2048" : input)

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
    Lx = parse(Float64, input == "" ? "100" : input)

    print("Ni (=10): ")
    input = readline()
    Ni = parse(Int64, input == "" ? "10" : input)

    surf = interface_prompt()

    println("Material above the interface (+)")
    above = material_prompt(default="vacuum")
    println("Material below the interface (-)")
    below = material_prompt(default="isotropic")
    
    print("Seed [Int64 >= 0] (= random seed): ")
    input = readline()
    seed = parse(Int64, input == "" ? "-1" : input)

    params = ParametersConfig(
        lambda=lambda,
        Nx=Nx,
        θs=angles,
        Lx=Lx,
        Ni=Ni,
        surf=surf,
        above=above,
        below=below,
        seed=seed,
    )

    default_name = splitext(basename(PROGRAM_FILE))[1]
    default_filename = "$(default_name)"
    
    print("Save as filename [=\"$(default_filename)\"]: ")
    input = readline()
    filename = input == "" ? default_filename : input
    
    save_parameters_config(path / filename, params)

    return nothing
end

function cli_help()
    println("""

    MDRC Commands:
    help - print this help
    new - create new run configuration
    input [path] - load configuration file or directory at \'path\', defaults to \'./input\'
    info [Int|name] - show information about configuration file
    run [Int|name] - run a configuration file
        out [label] - output file label following the run command
        iters [Int] - number of ensemble iters, default 100
        solver [full|reduced|hybrid] - solver method
    list - show list of loadable configurations in \'path\'
        and optionally run or show information

    """)
    return
end

function cli_list(path)
    files = get_files(path)
    for (idx, file) in enumerate(files)
        println("$idx: $file")
    end
    return
end

function cli_info(filepath)
    try
        params_config = load_parameters_config(filepath)
        display(params_config)
    catch
        try
            data = load_solver_data(filepath)
            display(data.params)
        catch
            error("File '$(filepath)' not found or file not valid")
        end
    end
    return
end

function cli_run(filepath, iters = 100, solver_type = :reduced)
    params_config = load_parameters_config(filepath)
    params = Parameters(params_config)
    display(params)

    @info "Initializing SolverData..."
    data = SolverData(params, iters, solver_type)
    @debug "CLI Init: SolverData complete"

    @info "Solving system of equations..."
    @time solve_ensemble!(data)

    if (idx = args_findoption("out")) !== nothing
        @assert length(ARGS) > idx "Missing argument after 'out'"
        fname = timestamp_suffix(ARGS[idx+1])
    else
        print("Custom label (datetime is appended): ")
        input = readline()
        fname = timestamp_suffix(input)
    end

    save_solver_data(DEFAULT_OUTPUT / fname, data)
    
    @info "Config and results saved to '$(DEFAULT_OUTPUT / fname)'"
    return
end

function cli_file(arg, argval, confpath, iters = 100, solver_type = :reduced)
    file_name = ""
    try
        files = get_files(confpath)
        file_name = files[parse(Int64, argval)]
    catch
        file_name = argval
    end

    filepath = confpath / file_name

    if (arg == "info")
        cli_info(filepath)
        exit(0)
    elseif (arg == "run")
        cli_run(filepath, iters, solver_type)
        exit(0)
    else
        error("Unknown single file argument: $(arg)")
        exit(0)
    end
end



function cli_directory(input_path, iters = 100, solver_type = :reduced)
    if (idx = args_findoption("list")) !== nothing
        cli_list(input_path)
        println("""
        Actions:
            info [Int] - show information about configuration file
            run [Int] - run configuration""")
        input = readline()

        arg, argval = split(input, " ")
        cli_file(arg, argval, input_path, iters, solver_type)
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
        input_path = ARGS[idx+1]
    end

    iters = 100
    if (idx = args_findoption("iters")) !== nothing
        @assert length(ARGS) > idx "Missing argument after 'iters'"
        iters = parse(Int64, ARGS[idx+1])
    end

    solver_type = :reduced
    if (idx = args_findoption("solver")) !== nothing
        @assert length(ARGS) > idx "Missing argument after 'solver'"
        solver_type = Symbol(ARGS[idx+1])
    end

    @debug input_path
    @debug isdir(input_path)
    @debug isfile(input_path)
    @debug pwd()

    if !isdir(input_path) && !isfile(input_path)
        @debug "Attempting relative path"
        input_path = DEFAULT_PATH / input_path
    end

    if (idx = args_findoption("info") !== nothing)
        if isfile(input_path)
            cli_info(input_path)
        elseif isdir(input_path)
            @assert length(ARGS) > idx "Missing argument after 'info'"
            cli_file("info", ARGS[idx+1], input_path, iters, solver_type)
        end
    elseif (idx = args_findoption("run") !== nothing)
        if isfile(input_path)
            cli_run(input_path, iters, solver_type)
        elseif isdir(input_path)
            @assert length(ARGS) > idx "Missing argument after 'run'"
            cli_file("run", ARGS[idx+1], input_path, iters, solver_type)
        end
    elseif isdir(input_path)
        cli_directory(input_path, iters, solver_type)
    end
end