#!/usr/bin/env julia

# Script to load a JLD2 file and dump it to a JSON file
using RayleighSolver

# Path to the source JLD2 file
input_file = "input/silver0.jld2"
# Path for the output JSON file
output_file = "input/silver0.json"

println("Loading parameters from $(input_file)...")
params = load_parameters(input_file)
println("Parameters loaded successfully.")

println("Saving parameters to $(output_file)...")
save_parameters_json(output_file, params)
println("Parameters saved to JSON successfully.")
