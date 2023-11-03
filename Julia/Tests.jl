include("SourceMatrices.jl")
using .SourceMatrices

using Plots
using BenchmarkTools

using SparseArrays
using LinearAlgebra
using StaticArrays
using LoopVectorization

function test_surfaces()
    # Units of µm here:
    L = 10.0
    δ = 30e-3
    a = 100e-3
    N = 1000
    M = 100

    @time se = SurfaceEnsemble(L, δ, a, N, M)

    display(se.surfs[1].ζ)

    plt = plot(
        size=(1200, 800),
        se.xs,
        se.surfs[1].ζ,
        ylims=(-0.1, 0.1))

    plot!(
        se.xs,
        se.surfs[M].ζ)

    display(plt)

    s = mean_slope(se.surfs[1])
    display("Numerical mean slope   $(s)")
    display("Analytical mean slope  $(√2*δ/a)")
end

# test_surfaces()

function test_matrix_gen()
    # Units of µm here:
    L = 10.0
    δ = 30e-3
    a = 100e-3
    N = 1000
    M = 100

    @time se = SurfaceEnsemble(L, δ, a, N, M)

    p = Params(1.0, 1e6, se.Δx)

    Asz = Matrix{ComplexF64}(undef, N, N)
    B = Matrix{ComplexF64}(undef, N, N)

    # Ast = spzeros(ComplexF64, N, N)


    bA1 = @benchmark create_A!($Asz, $se.surfs[1], $se.xs, $p) # 48 ms

    # bA2 = @benchmark create_A!($Ast, $se.surfs[1], $se.xs, $p) # 66 ms
    # bA3 = @benchmark create_A($se.surfs[1], $se.xs, $p) # 97 ms, 102 ms, 100 ms



    bB = @benchmark create_B!($B, $se.surfs[1], $se.xs, $p)

    create_A!(Asz, se.surfs[1], se.xs, p)
    create_B!(B, se.surfs[1], se.xs, p)

    display(Asz)
    display(B)



    # @code_warntype create_A!(A, se.surfs[1], se.xs, p)
    # @code_warntype create_B!(B, se.surfs[1], se.xs, p)

    display(bA1)
    # display(bA2)
    # display(bA3)

    display(bB)

end


test_surfaces()
