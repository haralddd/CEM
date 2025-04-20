"""
    function solve_single_full!(alloc::Preallocated, pre::Precomputed, data::SolverData)::Nothing

Calculates the preallocated surface integral.

# Arguments:
- `alloc`: [`Preallocated`](@ref) - Preallocated steps
- `pre`: [`Precomputed`](@ref) - Precomputed matrix elements
- `data`: [`SolverData`](@ref) - Contains the parameters, preallocated steps, and output structures

# Returns:
- Nothing
"""
function solve_single_full!(alloc::Preallocated, pre::Precomputed, data::SolverData)::Nothing
    @error "Not implemented for $(typeof(data))"
end
