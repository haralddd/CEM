# For MDRC plots
include("plot.jl")

const tol = 1.0e-10

function /(x::AbstractString, y::AbstractString)
    return joinpath(x, y)
end

const DEFAULT_PATH = joinpath(@__DIR__, "..")
const DEFAULT_OUTPUT = DEFAULT_PATH / "output"
const DEFAULT_PLOTDIR = DEFAULT_PATH / "plots"
mkpath(DEFAULT_PLOTDIR)

"""
    plot_lefty_data(filename::AbstractString, output_dir::AbstractString=DEFAULT_PLOTDIR)

Main function to generate specialized plots for left-handed materials.

# Arguments:
- `filename`: Name of the JLD2 file containing the simulation data
- `output_dir`: Directory to save plots to
"""
function plot_lefty_data(filename::AbstractString, output_dir::AbstractString=DEFAULT_PLOTDIR)
    solver_data = load_solver_data(filename)
    
    # Create output directory
    folder = output_dir / filename
    mkpath(folder)
    
    # Calculate MDRC and MDTC data
    mdrc_data = calc_mdrc(solver_data)
    θ0s = mdrc_data.θ0s
    θss = mdrc_data.θs
    
    # tuples of (ylabel, angle) 
    labels = [
        (L"\langle \mathrm{RC}\rangle_p", "rcp")
    ]

    fig_tcp = Figure(fontsize=24, size=(800, 600))
    fig_tcs = Figure(fontsize=24, size=(800, 600))



    for i in axes(cs, 2)
        
    end
    
    
    @info "Saved plots to $(folder)"
end

# Main execution when script is run directly
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        error("Usage: julia plot_lefty.jl <filename> [output_dir]")
    end
    
    filename = ARGS[1]
    output_dir = length(ARGS) > 1 ? ARGS[2] : DEFAULT_PLOTDIR
    
    plot_lefty_data(filename, output_dir)
end