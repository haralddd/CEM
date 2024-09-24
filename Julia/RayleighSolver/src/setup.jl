"c_0\\ 299_792_458.0 \\frac{m}{s}`` - Speed of light in a vacuum"
const c0 = 299_792_458.0
const FFT_Plan_t = FFTW.cFFTWPlan{ComplexF64,-1,true,1,Tuple{Int64}}

"""
Polarization of the primary field of the wave
"""
@enum Polarization p s

"""
    struct RayleighParams

Container of parameters which do not change during different realisations of the surface.
All lengths are scaled to ``x\\cdot\\frac{\\omega}{c_0}``, making length and wavenumber dimensionless.

# Fields
- `FT::FFTW.rFFTWPlan{T,K,inplace,dims,flags}`: Planned Fourier transform of surface points
- `xs::Vector{Float64}`: Surface points. Dimensionless, in units of ``x\\cdot\\frac{\\omega}{c_0}``
- `ps::Vector{Float64}`: Discretized transmitted wave numbers. Dimensionless, in units of ``k\\cdot\\frac{c_0}{\\omega}``
- `qs::Vector{Float64}`: Discretized reflected wave numbers. Dimensionless, in units of ``k\\cdot\\frac{c_0}{\\omega}``
- `ks::Vector{Float64}`: Discretized incident wave numbers. Dimensionless, in units of ``k\\cdot\\frac{c_0}{\\omega}``
- `kis::Vector{Int64}`: Indices of k in qs
- `nu::Polarization`: [`Polarization`](@ref) - Enum of [p, s] type primary field
- `eps::ComplexF64`: Permittivity of the scattering medium
- `mu::ComplexF64`: Permeability of the scattering medium
- `lambda::Float64`: Wavelength in [m]
- `omega::Float64`: Angular frequency
- `Q::Int64`: Truncated wave number multiplum, ``q = (-\\infty,\\infty) \\rightarrow q_n \\in (-\\frac{Q}{2}, \\frac{Q}{2})``
- `dq::Float64`: Wave number spacing [unitless]
- `Lx::Float64`: Surface length [unitless]
- `dx::Float64`: Surface point spacing [unitless]
- `Ni::Int`: Number of terms in the Taylor expansion of ``I(\\gamma | q)``
- `Nq::Int`: Number of discretised reflected wave numbers
- `Nx::Int`: Number of surface points
- `surf::SurfaceParams`: [SurfaceParams](@ref) - Surface parameter pack
- `seed::Int`: Random seed for the algorithm
- `rng::Xoshiro`: Random number generator with the given seed
"""
struct RayleighParams{SurfType<:SurfaceParams}
    FFT::FFT_Plan_t

    xs::Vector{Float64}
    xks::Vector{Float64}
    
    ps::Vector{Float64}
    qs::Vector{Float64}
    ks::Vector{Float64}
    kis::Vector{Int}

    nu::Polarization
    eps::ComplexF64
    mu::ComplexF64
    lambda::Float64
    omega::Float64

    Q::Int64
    dq::Float64
    Lx::Float64
    dx::Float64

    Ni::Int
    Nq::Int
    Nx::Int

    surf::SurfType
    seed::Int
    rng::Xoshiro
    function RayleighParams(;
        nu::Polarization=p,
        eps=2.25,
        mu=1.0,
        lambda=600e-9,
        Q=4,
        Nq=127,
        ks=[sind(10.0)],
        L=10.0e-6,
        Ni=10,
        surf::SurfType,
        seed=-1,
        rescale=true
    )::RayleighParams where SurfType<:SurfaceParams
    
        K = 2π / lambda
        omega = c0 * K
    
        if imag(eps * mu) == 0.0
            println("epsmu is purely real, adding a 1e-4i to eps to avoid singularities")
            eps += 1e-4im # Add a small imaginary part to avoid singularities
        end
    
    
        # Assertions and warnings
        @assert L / lambda > 10.0 "Surface length must be much larger than the wavelength, L >> lambda, but is L:$L and lambda:$lambda."
        @assert Q > 2 "Q must be greater than 2, but is $Q. 4 is recommended."
        @assert Nq > 2 "Nq must be greater than 2, but is $Nq."
    
        if rescale
            Lx = L * omega / c0
            new_surf = scale(surf, omega/c0)
        else
            Lx = L
            new_surf = surf
        end
        Nx = 2 * Nq
        dq = Q / (Nq-1)
    
        dx = Lx / (Nx-1)
        xs = -Lx/2:dx:Lx/2
        xks = fftfreq(Nx, 2pi / dx)
    
        ps = -Q/2:dq:Q/2
        qs = -Q/2:dq:Q/2
    
        kis = [searchsortedfirst(qs, k) for k in ks] |> collect
        ks = qs[kis]
    
        seed = seed == -1 ? rand(0:typemax(Int)) : seed
        new{SurfType}(plan_fft!(similar(xs, ComplexF64)),
            xs, xks, ps, qs, ks, kis,
            nu, eps, mu, lambda, omega,
            Q, dq, Lx, dx,
            Ni, Nq, Nx, new_surf, seed, Xoshiro(seed))
    end
    
    function RayleighParams(serial_file::String)::RayleighParams
        open(serial_file, "r") do io
            obj = deserialize(io)
            return obj
        end
    end

end



"""
    get_angles(rp::RayleighParams)::Vector{Float64}

Returns the incident angles of `rp` in degrees.
"""
function get_angles(rp::RayleighParams)::Vector{Float64}
    return asind.(rp.qs[rp.kis])
end

function get_scale(rp::RayleighParams)::Float64
    return rp.omega / c0
end

function scaled_params(rp::RayleighParams)::Dict
    s = get_scale(rp)
    return Dict(
        :nu => rp.nu,
        :eps => rp.eps,
        :mu => rp.mu,
        :lambda => rp.lambda,
        :omega => rp.omega,
        :Q => rp.Q,
        :ks => rp.ks * s,
        :Nq => rp.Nq,
        :Lx => rp.Lx / s,
        :Ni => rp.Ni,
        :surf => scale(rp.surf, 1/s),
    )
end