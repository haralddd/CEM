"""
Functions for T-matrix calculations in the CEM framework.
"""

using FFTW
using LinearAlgebra

function In(ζ, γ, q, n, fft)
    return (-1im * γ)^n / factorial(n)
end

function compute_T(R, data)
    return R .* sum()    
end
