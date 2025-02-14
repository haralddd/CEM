using CEM

# Fixed parameters for all configurations
const λ = 589.3e-9  # nm, sodium D-line
const Nx = 2048
const Lx = 100 * λ
const Ni = 10

# Parameter ranges to sweep
angles = 0.0:10.0:60.0  # degrees
heights = [50e-9, 100e-9, 150e-9]  # different surface heights
seeds = 1:5  # different random seeds for each configuration

# Create output directory if it doesn't exist
output_dir = joinpath("input", "sweep")
mkpath(output_dir)

# Generate configurations for all combinations
for h in heights, θ in angles, seed in seeds
    # Create surface with given height
    surf = GaussianSurface(h, Lx/20)  # correlation length = Lx/20
    
    # Create materials (example with air/glass interface)
    above = Material(1.0 + 0.0im)  # air
    below = Material(1.5 + 0.0im)  # glass
    
    # Create parameters
    params = Parameters(
        lambda=λ,
        Nx=Nx,
        θs=[θ],  # single angle per configuration
        Lx=Lx,
        Ni=Ni,
        surf=surf,
        above=above,
        below=below,
        seed=seed,
        rescale=true,
    )
    
    # Generate filename based on parameters
    filename = "config_h$(Int(round(h*1e9)))nm_theta$(Int(θ))deg_seed$(seed).jld2"
    
    # Save configuration
    save_spa_config(joinpath(output_dir, filename), params)
    println("Saved configuration: $filename")
end
