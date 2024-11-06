"c_0\\ 299_792_458.0 \\frac{m}{s}`` - Speed of light in a vacuum"
const c0 = 299_792_458.0
const FFT_Plan_t = FFTW.cFFTWPlan{ComplexF64,-1,true,1,Tuple{Int64}}
const IFFT_Plan_t = AbstractFFTs.ScaledPlan{ComplexF64,FFTW.cFFTWPlan{ComplexF64,1,true,1,UnitRange{Int64}},Float64}


"""
    SimParams{SurfT<:RandomSurface, Above<:Material, Below<:Material}(

    )

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
- `surf::RandomSurface`: [RandomSurface](@ref) - Surface parameter pack
- `seed::Int`: Random seed for the algorithm
- `rng::Xoshiro`: Random number generator with the given seed
"""
struct SimParams{SurfT<:RandomSurface, AboveT<:Material, BelowT<:Material}
    FFT::FFT_Plan_t
    IFFT::IFFT_Plan_t

    xs::Vector{Float64}
    xks::Vector{Float64}

    ps::Vector{Float64}
    qs::Vector{Float64}
    ks::Vector{Float64}
    kis::Vector{Int}
    rev_qis::Vector{Int}
    rev_kis::Vector{Int}

    above::AboveT
    below::BelowT
    lambda::Float64
    omega::Float64

    Q::Int64
    dq::Float64
    Lx::Float64
    dx::Float64

    Ni::Int
    Nq::Int
    Nx::Int

    surf::SurfT
    seed::Int
    rng::Xoshiro
end

function SimParams{ST,AT,BT}(;
    lambda=600e-9,
    Q=4,
    Nq=127,
    ks=[sind(10.0)],
    Lx=10.0e-6,
    Ni=10,
    surf::ST = FlatSurface(),
    above::AT = Vacuum(),
    below::BT = Isotropic(2.25 + 1e-6im, 1.0),
    seed=-1,
    rescale=true
) where {
    ST<:RandomSurface,
    AT<:Material,
    BT<:Material}

    K = 2π / lambda
    omega = c0 * K

    # Assertions and warnings
    @assert Lx / lambda > 10.0 "Surface length must be much larger than the wavelength, L >> lambda, but is L:$L and lambda:$lambda."
    @assert Q > 2 "Q must be greater than 2, but is $Q. 4 is recommended."
    @assert Nq > 2 "Nq must be greater than 2, but is $Nq."

    if typeof(below) == Isotropic
        if imag(below.eps) + imag(below.mu) ≈ 0
            @warn "Material below has no loss, adding small imaginary component to avoid singularities."
            below = Isotropic(below.eps + 1e-4im, below.mu)
        end
    end

    if rescale
        Lx = Lx * omega / c0
        new_surf = scale(surf, omega / c0)
    else
        Lx = Lx
        new_surf = surf
    end

    Nx = 2 * Nq
    dq = Q / (Nq - 1)

    dx = Lx / (Nx - 1)
    xs = -Lx/2:dx:Lx/2
    xks = fftfreq(Nx, 2pi / dx)

    ps = -Q/2:dq:Q/2
    qs = -Q/2:dq:Q/2

    kis = [searchsortedfirst(qs, k) for k in ks] |> collect
    ks = qs[kis]

    # To access the the Fourier transform of the surface integral I(γ, q)
    # we must access the pattern ζ(p-q), so make a reverse index range for q
    # This works only if qs is symmetric (-Q/2:0:Q/2)
    rev_qis = reverse(eachindex(qs))
    @assert all(qs[rev_qis] .== .-qs)

    # for k we must find the index of the corresponding q with opposite sign
    rev_kis = [searchsortedlast(qs, -k) for k in ks] |> collect
    @assert all(qs[rev_kis] .== -qs[kis])

    @debug "SimParams"
    @debug Lx
    @debug new_surf
    @debug xks
    @debug xs
    @debug ps
    @debug qs
    @debug ks
    @debug kis
    @debug rev_qis
    @debug rev_kis

    seed = seed < 0 ? rand(0:typemax(Int)) : seed
    FFTplan = plan_fft!(similar(xs, ComplexF64))
    IFFTplan = plan_ifft!(similar(xs, ComplexF64))
    @debug typeof(FFTplan)
    @debug typeof(IFFTplan)
    SimParams{ST,AT,BT}(FFTplan, IFFTplan,
        xs, xks, ps, qs, ks, kis, rev_qis, rev_kis,
        above, below, lambda, omega,
        Q, dq, Lx, dx,
        Ni, Nq, Nx, new_surf, seed, Xoshiro(seed))
end

function SimParams(;
    lambda=600e-9,
    Q=4,
    Nq=127,
    ks=[sind(10.0)],
    Lx=10.0e-6,
    Ni=10,
    surf = FlatSurface(),
    above = Vacuum(),
    below = Isotropic(2.25 + 1e-6im, 1.0),
    seed=-1,
    rescale=true
)::SimParams
    return SimParams{typeof(surf), typeof(above), typeof(below)}(
        lambda=lambda,
        Q=Q,
        Nq=Nq,
        ks=ks,
        Lx=Lx,
        Ni=Ni,
        surf=surf,
        above=above,
        below=below,
        seed=seed,
        rescale=rescale
    )
end