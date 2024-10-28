# using MKL
push!(LOAD_PATH, "$(@__DIR__)/RayleighSolver/")
using RayleighSolver
using Dates
using Statistics

function /(x::String, y::String)
    return joinpath(x, y)
end

const DEFAULT_PATH = "$(@__DIR__)/"
const DEFAULT_INPUT = DEFAULT_PATH / "input"
const DEFAULT_OUTPUT = DEFAULT_PATH / "output"

function timestamp_suffix(str="out")
    return str * "_" * Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")
end

function args_findoption(option)::Union{Nothing,Int}
    return findfirst(isequal(option), ARGS)
end

function get_jld_files(path)
    return filter(x -> endswith(x, ".jld2"), readdir(path))
end


function cmd_help()
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

function cmd_list(path)
    files = get_files(path)
    for (idx, file) in enumerate(files)
        println("$idx: $file")
    end
    return
end

function cmd_info(filepath)
    try
        spa = load_spa_config(filepath)
        iters = load_ensemble_iters(filepath)
        display(spa)
        display("Ensemble iterations: $iters")
    catch
        error("File '$(filepath)' not found or not valid")
        exit(0)
    end
    return
end

function cmd_run(filepath)
    cmd_info(filepath)
    spa = load_spa_config(filepath)
    iters = load_ensemble_iters(filepath)
    prealloc = SimPrealloc(spa)
    qs, coh, incoh = solve_MDRC!(spa, prealloc, iters)
    return
end

function cmd_single_file(arg, argval, confpath)
    file_name = ""
    try
        files = get_files(confpath)
        file_name = files[parse(Int64, argval)]
    catch
        file_name = argval
    end

    filepath = confpath / file_name

    if (arg == "info")
        cmd_info(filepath)
        exit(0)
    elseif (arg == "run")
        cmd_run(filepath)
        exit(0)
    else
        error("Unknown single file argument: $(arg)")
        exit(0)
    end
end



function cmd_directory(input_path)
    if (idx = args_hasoption("list")) !== nothing
        cmd_list(input_path)
        println("""
        Actions:
            info [Int] - show information about configuration file
            run [Int] - run configuration""")
        input = readline()

        arg, argval = split(input, " ")
        cmd_file(arg, argval, input_path)
    else
        error("Unknown action following input path: $(input_path)")
        exit(0)
    end
end



function cmd_main()
    if length(ARGS) < 1 || ARGS[1] == "help"
        cmd_help()
        exit(0)
    end

    if (idx = args_hasoption("new")) !== nothing
        config_creation_prompt()
        exit(0)
    end

    input_path = "input"

    if (idx = args_hasoption("input")) !== nothing
        @assert length(ARGS) > idx "Missing argument after 'input'"
        try
            input_path = parse(String, ARGS[idx+1])
        catch
            error("Invalid positional value after 'input': $(ARGS[idx+1])")
            exit(0)
        end
    end

    if (idx = args_hasoption("info") !== nothing 
        || (idx = args_hasoption("run")) !== nothing)
        @assert length(ARGS) > idx "Missing argument after 'info' or 'run'"
        if isfile(input_path)
            cmd_file(ARGS[idx], ARGS[idx+1], input_path)
        end
    elseif isdir(input_path)
        cmd_directory(input_path)
    end
end

if (abspath(PROGRAM_FILE) == @__FILE__) 
    cmd_main()
end