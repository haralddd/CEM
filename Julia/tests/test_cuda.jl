using CUDA
using CUDA.CUFFT
using FFTW
using Random
using LinearAlgebra

const N = 2048

function test_cuda()
    CUDA.allowscalar(false)
    dtype = ComplexF32
    # dtype = Float64
    CUDA.versioninfo()

    A = CUDA.randn(dtype, (N, N))
    B = CUDA.randn(dtype, (N, N))
    C = similar(B)

    b = CUDA.randn(dtype, (N))
    x = CUDA.randn(dtype, (N))
    F = CUDA.CUFFT.plan_fft(b)


    @info "CUDA FFT and Ax=b"
    CUDA.@time begin
        for i in 1:100
            CUDA.randn!(x)
            b = F * x

            C .= A .* B
            x = A \ b
        end
        nothing
    end
    nothing
end

function test_cpu()
    dtype = ComplexF32
    # dtype = Float64
    A = ones((N,N)) .+ Random.randn(dtype, (N, N))
    B = ones((N,N)) .+ Random.randn(dtype, (N, N))
    C = similar(B)

    b = Random.randn(dtype, (N))
    x = similar(b)
    F = FFTW.plan_fft(b)


    @info "CPU FFT and Ax=b"
    @time begin
        for i in 1:100
            Random.randn!(x)
            b = F * x

            C .= A .* B
            x = A \ b
        end
        nothing
    end
    nothing
    display(typeof(x))
end

test_cuda()
test_cpu()