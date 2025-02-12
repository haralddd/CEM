module DiscreteDiff
export d1o2, d2o2

function d1o2(x::VecType, Δx::T)::VecType where {T<:Real,VecType<:AbstractVector{T}}
    # https://en.wikipedia.org/wiki/Finite_difference_coefficient
    dx = similar(x)
    N = length(dx)
    @assert N > 2

    # Forward on start point to order O(h²)
    dx[1] = -1.5 * x[1] + 2.0 * x[2] - 0.5 * x[3]
    # Backward on end points to order O(h²)
    dx[end] = 1.5 * x[end] - 2.0 * x[end-1] + 0.5 * x[end-2]

    # Central on the rest to order O(h²)
    for i in 2:N-1
        dx[i] = -0.5 * x[i-1] + 0.5 * x[i+1]
    end

    return dx ./ Δx
end
function d2o2(x::VecType, Δx::T)::VecType where {T<:Real,VecType<:AbstractVector{T}}
    # https://en.wikipedia.org/wiki/Finite_difference_coefficient
    dx2 = similar(x)
    N = length(dx2)
    @assert N > 3

    # Forward on start point to order O(h²)
    dx2[1] = 2.0 * x[1] - 5.0 * x[2] + 4.0 * x[3] - 1.0 * x[4]
    # Backward on end points to order O(h²)
    dx2[end] = 2.0 * x[end] - 5.0 * x[end-1] + 4.0 * x[end-2] - 1.0 * x[end-3]

    # Central on the rest to order O(h²)
    for i in 2:N-1
        dx2[i] = 1.0 * x[i-1] - 2.0 * x[i] + 1.0 * x[i+1]
    end

    return dx2 ./ Δx^2
end
end # module