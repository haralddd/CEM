# using MKL
push!(LOAD_PATH, "Julia/RayleighSolver/")
using RayleighSolver
using Dates

function setup_dir(str="data")
    display("Adding suffix to dir")
    str *= "_" * Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")
    mkdir(str)
    display("New dir: \"$str\"")
    return str
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("Running as script...")
    print("Make new run configuration [y] or load existing [n]? ")
    if readline() == "y"
        rp, sp, generator! = make_solver_config()
    elseif readline() == "n"
        print("Enter config file name (./input/{filename}): ")
        rp, sp, generator! = load_solver_config(readline())
    end
else
    println("Running as module...")
    rp, sp, generator! = make_solver_config()
    exit(1)
end

print("Save data in [./data/{simulation_name}]: ")
run_dir = "data/" * readline()
run_dir = setup_dir(run_dir)
save_solver_config(run_dir * "/config.txt", rp)

coh, incoh = run_threaded(rp; N_ens=N, θ=θ)

# Save data
open(run_dir * "/coh.bin", "w") do io
    write(io, coh)
end
open(run_dir * "/incoh.bin", "w") do io
    write(io, incoh)
end