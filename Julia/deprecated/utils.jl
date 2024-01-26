#=
    This file contains utility functions that may or may not be used in the
    main code. They are here for convenience.
=#

function is_power_two(n::Int)::Bool
    # Check if n is a power of 2
    (n & (n - 1)) == 0
end


function trapz(xs, ys)
    # Simple trapezoidal integration with variable step size
    @assert length(xs) == length(ys)
    res = 0.0

    for i in eachindex(xs)
        i == 1 && continue
        res += 0.5 * abs(xs[i] - xs[i-1]) * (ys[i] + ys[i-1])
    end
    return res
end