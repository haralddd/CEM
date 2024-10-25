"""
    SimPrealloc{N<:Integer}(Mpq::Matrix{ComplexF64}, Npk::Matrix{ComplexF64}, Fys::Vector{ComplexF64}, sFys::Vector{ComplexF64}, ys::Vector{Float64})

SimPrealloc{N<:Integer} is a struct which contains all
preallocated matrices which are used during
the calculation of the surface integral, i.e. the array values are mutated in place
Note that all members are uninitialized, and must be initialized after making the struct

# Fields

- `Mpq::Matrix{ComplexF64}`: Matrix of the Mpq coefficients (A)
- `Npk::Matrix{ComplexF64}`: Vector (for all k) of Vectors of the Npk coefficients (b)
- `Fys::Vector{ComplexF64}`: Fourier transform of surface heights, prealloc
- `sFys::Vector{ComplexF64}`: Shifted Fourier transform of surface heights, prealloc
- `Z::Vector{ComplexF64}`: Preallocated computation step in surface generation
- `ys::Vector{Float64}`: Surface heights
"""
mutable struct SimPrealloc
    Mpq::Matrix{ComplexF64}
    Npk::Matrix{ComplexF64}
    FM::LU{ComplexF64, Matrix{ComplexF64}, Vector{Int64}}
    Fys::Vector{ComplexF64}
    sFys::Vector{ComplexF64}
    Z::Vector{ComplexF64}
    ys::Vector{Float64}

    function SimPrealloc(Nq::Integer, Nk::Integer)
        ys = Vector{Float64}(undef, 2Nq)
        Fys = similar(ys, ComplexF64)
        sFys = similar(Fys)
        Mpq = Matrix{ComplexF64}(undef, Nq, Nq)
        Npk = Matrix{ComplexF64}(undef, Nq, Nk)
        FM = lu(Mpq; check=false)
        Z = similar(Fys)

        new(Mpq, Npk, FM, Fys, sFys, Z, ys)
    end
    function SimPrealloc(spa::SimParams)
        SimPrealloc(spa.Nq, length(spa.ks))
    end
end