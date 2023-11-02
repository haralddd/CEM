#=
    This file contains utility functions that may or may not be used in the
    main code. They are here for convenience.
=#

function is_power_two(n::Int)::Bool
    # Check if n is a power of 2
    (n & (n - 1)) == 0
end