include("test_surface.jl")
include("test_fresnel.jl")
include("test_unitary.jl")
include("test_solver.jl")

# Run surface tests
test_gaussian()
test_rect()

# Run fresnel tests
test_fresnel()

# Run unitary tests
test_unitary()

# Run complete tests
test_solver()