# # Run surface tests
# begin
#     include("test_surface.jl")
#     test_gaussian()
#     test_gaussian2()
#     test_gaussian3()
# end
# test_rect()

# # Run fresnel tests
# include("test_fresnel.jl")
# test_fresnel()

# # Run unitary tests
# include("test_unitary.jl")
# test_unitary()

# Run complete solver tests
# include("test_solver.jl")

# test_reciprocity()
# test_symmetry_isotropic()
# test_solver() 
# profile_solver_components()
# test_crystal_precompute()

include("benchmarks.jl")
profile_gaussian_surfacegen()
profile_rectangular_surfacegen()
profile_isotropic_precompute()
profile_isotropic_single_solve()
