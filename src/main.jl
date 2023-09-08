include("SourceMatrices.jl")
using .SourceMatrices

using QuadGK

struct BeamParams
    w::Float64 # Half-width
    ω::Float64 # Frequency region width
end

function α0(q, ω)::ComplexF64

    return abs(q) < ω / c ?
           √((ω / c)^2 - q^2) :
           1im * √(q^2 - (ω / c)^2)

end

function Fg(k::Float64, p::BeamParams)
    # Gaussian finite beam envelope

    w = p.w
    ω = p.ω

    return w * ω / (2 * √π * c * α0(k, ω)) * exp(-w^2 * ω^2 / (4 * c^2) * (asin(k * c / ω) - θ)^2)

end

function Φ_inc_g(x1, x3, ω, p::BeamParams)
    #=
    Computes finite width incident wave with Gaussian envelope

    =#

    # Kernel of integration of the incident wavepacket
    ker(k) = Fg(k, p) * exp(1im * k * x1 - 1im * α0(k, ω) * x3) / 2π

    # Endpoints: 
    # TODO: Check complex integration(up to but not including ±ω/c, due to div by zero in F(k))
    #   If numerical, we may just define a wave Φ directly without this analytic expression?
    a = -ω / c
    b = ω / c

    # Quadrature integration of the source function
    return quadgk(ker, a, b)
end



function solver(s::Surface, ξ::Vector{Float64}, p⁺::Params, p⁻::Params)
    N = length(ξ)

    # Source matrices for outer medium
    A⁺ = Matrix{ComplexF64}(undef, N, N)
    B⁺ = Matrix{ComplexF64}(undef, N, N)
    create_A!(A⁺, s, ξ, p⁺)
    create_B!(B⁺, s, ξ, p⁺)

    # Source matrices for inner medium
    A⁻ = Matrix{ComplexF64}(undef, N, N)
    B⁻ = Matrix{ComplexF64}(undef, N, N)
    create_A!(A⁻, s, ξ, p⁻)
    create_B!(B⁻, s, ξ, p⁻)


    Φ = Matrix{ComplexF64}(undef, N, N)
    mask = Matrix{Bool}(undef, size(Φ))
    for I in eachindex(IndexCartesian(), mask)
        i, j = I.I
        mask[I] = s.ζ[i] > ξ[j]
    end
    Φ⁺ = Φ[mask]
    Φ⁻ = Φ[.!mask]

    display(mask)
    display(typeof(Φ⁺))
    display(typeof(Φ⁻))




    # Solve the source function matrix equations
    # F_new = Vector{N,ComplexF64}(undef)
    # F_inc = SizedVector{N,ComplexF64}(undef)


    nothing
end

function main()
    # Units of µm here:
    L = 10.0
    δ = 30e-3
    a = 100e-3
    N = 1000
    M = 100

    se = SurfaceEnsemble(L, δ, a, N, M)


    p⁺ = Params(
        ε=ε_AIR,
        ω=1e6,
        Δξ=se.Δx)


    p⁻ = Params(
        ε=1.0,
        ω=1e6,
        Δξ=se.Δx)


    solver(se.surfs[1], se.xs, p⁺, p⁻)
end

main()
