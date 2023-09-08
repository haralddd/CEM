module SourceMatrices
export create_A!, create_B!, SourceParams

# TODO: Create specialized symmetric matrix with static size (such as SizedMatrix from StaticArrays.jl)

struct SourceParams
    # Parameter pack for source matrix generation
    ε::Float64
    ω::Float64
    Δξ::Float64
end

SourceParams(; ε, ω, Δξ) = SourceParams(ε, ω, Δξ)


function create_A!(A::AbstractMatrix{ComplexF64}, s::Surface, ξ::AbstractVector{Float64}, p::SourceParams)::Nothing
    c = 299_792_458.0 # speed of light m/s
    γ = 0.5772156649015328606065 # Euler constant

    ε = p.ε
    ω = p.ω
    Δξ = p.Δξ
    N, M = size(A)
    @assert M == N "Expected square matrix"

    γf(m::Int)::Float64 =
        √(1.0 + s.ζ_dot[m]^2)

    # Hankell function of first kind approximation around small z
    H₁(z::ComplexF64)::ComplexF64 =
        -2.0im / (π * z^2) + 1.0im / π * (log(0.5 * z) + γ + 0.5) - 0.5

    χ(m::Int, n::Int)::ComplexF64 =
        √ε * ω / c * √((ξ[m] - ξ[n])^2 + (s.ζ[m] - s.ζ[n])^2)

    # Kernel of A
    A_ker(m::Int, n::Int)::ComplexF64 = (
        -0.25im * ε * (ω / c)^2 * H₁(χ(m, n)) / χ(m, n) *
        ((ξ[m] - ξ[n]) * s.ζ_dot[n] - (s.ζ[m] - s.ζ_dot[n]))
    )

    # Matrix elements
    A_mn(m::Int, n::Int)::ComplexF64 =
        Δξ * A_ker(m, n)

    A_mm(m::Int)::ComplexF64 =
        0.5 + Δξ * s.ζ_ddot[m] / (4π * (γf(m))^2)

    for n in 1:N
        # Upper triangular elements
        for m in 1:n-1
            A[m, n] = A_mn(m, n)
        end

        A[n, n] = A_mm(n)

        # Lower triangular elements
        for m in n+1:N
            A[m, n] = A_mn(m, n)
        end
    end

    return nothing
end


function create_B!(B::AbstractMatrix{ComplexF64}, s::Surface, ξ::AbstractVector{Float64}, p::Params)::Nothing
    c = 299_792_458.0 # speed of light m/s
    γ = 0.5772156649015328606065 # Euler constant

    ε = p.ε
    ω = p.ω
    Δξ = p.Δξ

    N, M = size(B)
    @assert M == N "Expected square matrix"

    γf(m::Int)::Float64 =
        √(1.0 + s.ζ_dot[m]^2)

    # Hankell function of first kind approximations around small z
    H₀(z::Float64)::ComplexF64 =
        2.0im / π * (log(0.5 * z) + γ) + 1.0

    χ(m::Int, n::Int)::Float64 =
        √ε * ω / c * √((ξ[m] - ξ[n])^2 + (s.ζ[m] - s.ζ[n])^2)

    # Kernel of B
    B_ker(m::Int, n::Int)::ComplexF64 = (
        -0.25im * H₀(χ(m, n))
    )

    B_mn(m::Int, n::Int)::ComplexF64 =
        Δξ * B_ker(m, n)

    B_mm(m::Int)::ComplexF64 =
        -0.25im * Δξ * H₀(√ε * ω * γf(m) * Δξ / (2 * exp(1) * c))

    for n in 1:N
        # TODO: These elements are always symmetric, use Symmetric on sparse array. See StaticArrays.jl: SHermitianCompact and modify
        for m in 1:n-1
            # Upper triangular elements
            B[m, n] = B_mn(m, n)
        end

        B[n, n] = B_mm(n)

        for m in n+1:N
            # Upper triangular elements
            B[m, n] = B_mn(m, n)
        end
    end

    return nothing
end

end #module SourceMatrices