
using Pkg
Pkg.add("CUDA")
using CUDA
Pkg.test("CUDA")

function test()
    A = CuArray(randn(1024, 1024))
    b = CuArray(randn(1024))
    x = A / b

    return x
end