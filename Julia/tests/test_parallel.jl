
using CUDA
using LoopVectorization
using BenchmarkTools
using CUDA.CUFFT


function test_parallel()
    # CUDA.versioninfo()

    N = 1024
    A = CUDA.randn(ComplexF64, (N, N))
    b = CUDA.randn(ComplexF64, (N))
    x = similar(b)
    F = CUDA.CUFFT.plan_fft!(b, (1,))
    # IFFT = plan_ifft!(b)

    @btime $F*$b

    @btime $x = $A \ $b


end

test_parallel()