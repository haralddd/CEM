
using FFTW
using Random


function test_rfft(N = 10000)
    Nx = 2048
    ys = zeros(Nx)
    Fys = rfft(ys)

    RFFT = plan_rfft(ys)
    IRFFT = plan_irfft(Fys, Nx)

    rng = Xoshiro(1234)
    A = randn(length(Fys))

    @time for n in 1:N
        randn!(rng, ys)
        Fys .= RFFT * ys
        Fys .*= A
        ys .= IRFFT * Fys
    end
    return ys
end


function test_fft_inplace(N = 10000)
    Nx = 2048
    ys = zeros(Nx)
    Fys = zeros(ComplexF64, Nx)

    FFT = plan_fft!(Fys)
    IFFT = plan_ifft!(Fys)
    rng = Xoshiro(1234)
    A = randn(Nx)

    @time for n in 1:N
        randn!(rng, ys)
        Fys .= ys
        FFT * Fys
        Fys .*= A
        IFFT * Fys
        ys = real(Fys)
    end
    return ys
end

function test_fft_inplace_shifts(N=10000)
    Nx = 2048
    ys = zeros(Nx)
    Fys = zeros(ComplexF64, Nx)
    sFys = similar(Fys)

    FFT = plan_fft!(Fys)
    IFFT = plan_ifft!(Fys)
    rng = Xoshiro(1234)
    A = randn(Nx)

    @time for n in 1:N
        randn!(rng, ys)
        Fys .= ys

        FFT * Fys
        fftshift!(sFys, Fys)

        sFys .*= A

        ifftshift!(Fys, sFys)
        IFFT * Fys

        ys = real(Fys)
    end
    return ys
end

y1 = test_rfft()
y2 = test_fft_inplace()
y3 = test_fft_inplace_shifts()