"c_0\\ 299_792_458.0 \\frac{m}{s}`` - Speed of light in a vacuum"
const c0 = 299_792_458.0
const FFT_Plan_t = FFTW.cFFTWPlan{ComplexF64,-1,true,1,Tuple{Int64}}
const IFFT_Plan_t = AbstractFFTs.ScaledPlan{ComplexF64,FFTW.cFFTWPlan{ComplexF64,1,true,1,UnitRange{Int64}},Float64}


struct ParametersConfig{ST<:RandomSurface,AT<:Material,BT<:Material}
    lambda::Float64
    Lx::Float64
    Nx::Int
    θs::Vector{Float64}
    Ni::Int
    surf::ST
    above::AT
    below::BT
    seed::Int

    function ParametersConfig(;
        lambda=632.8e-9,
        Lx=100,
        Nx=2048,
        θs=[0.0, 10.0, 20.0],
        Ni=10,
        surf::ST=GaussianSurface(10.0e-9, 100.0e-9),
        above::AT=Vacuum(),
        below::BT=Isotropic(2.25 + 1e-6im, 1.0),
        seed=-1
    ) where {
        ST<:RandomSurface,
        AT<:Material,
        BT<:Material}
        return new{ST,AT,BT}(lambda, Lx, Nx, θs, Ni, surf, above, below, seed)
    end
end

"""
    Parameters{SurfT<:RandomSurface, Above<:Material, Below<:Material}(

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
struct Parameters{SurfT<:RandomSurface,AboveT<:Material,BelowT<:Material}
    lambda::Float64
    Lx::Float64
    Nx::Int
    Nq::Int
    Ni::Int
    dx::Float64
    dq::Float64

    xs::Vector{Float64}
    xks::Vector{Float64}
    ps::Vector{Float64}
    qs::Vector{Float64}

    θs::Vector{Float64}
    ks::Vector{Float64}
    kis::Vector{Int}
    sFys_pqidxs::Matrix{Int}

    surf::SurfT
    above::AboveT
    below::BelowT

    seed::Int
    rng::Xoshiro
    FFT::FFT_Plan_t
    IFFT::IFFT_Plan_t

    function Parameters(config::ParametersConfig{ST,AT,BT}) where {
        ST<:RandomSurface,
        AT<:Material,
        BT<:Material}

        lambda = config.lambda
        Lx = config.Lx
        Nx = config.Nx
        θs = config.θs
        Ni = config.Ni
        surf = config.surf
        above = config.above
        below = config.below
        seed = config.seed

        k0 = 2π / lambda # = ω/c, inverse length scale being used for dispersion relations and in real space

        # Assertions and warnings
        @assert Lx > 10.0 "Surface length must be much larger than the wavelength, Lx >> lambda, but is Lx:$Lx and lambda:$lambda."

        # Rescale parameters when constructed from the config
        Lx = Lx * k0 * lambda
        surf = scale(surf, k0)

        dx = Lx / Nx
        xs = (-Lx+dx)/2:dx:(Lx-dx)/2

        @debug "xs has zero: $(any(xs .== 0.0))"

        Nq = Nx ÷ 2
        Q = π / dx
        dq = 2π / Lx
        xks = fftfreq(Nx, 2pi / dx)

        @assert Nx % 2 == 0 "Nx should be divisible by 2 to make Nq."

        @debug "Q: $Q"
        @debug "xks has zero: $(any(xks .== 0.0))"
        @assert Q > 4 "Q (= $(Q)) should be greater than 4, which is determined by spatial resolution, Q=π/dx, with dx = Lx/Nx"

        if BT == Uniaxial
            no = sqrt(below.mu_para * below.eps_para)
            ne = sqrt(below.mu_perp * below.eps_perp)
            @debug "no: $no"
            @debug "ne: $ne"
            if Q / 2 < real(no) || Q / 2 < real(ne)
                @error "Q/2 ($Q/2) is smaller than no=$(no) or ne=$(ne). Not all transmission coefficients can be resolved."
            end
        end

        ps = -Q/2:dq:Q/2
        qs = -Q/2:dq:Q/2

        @debug "ps has zero: $(any(ps .== 0.0))"
        @debug "ps: $ps, qs: $qs"
        @debug "Nq: $Nq, len ps: $(length(ps)), len qs: $(length(qs)), Nx: $(length(xs))"
        @debug "dq: $dq, Δq: $(qs[2] - qs[1]) vs Δxks: $(xks[2] - xks[1])"
        @debug "min(xks): $(minimum(xks)), qs,ps[1]: $(qs[1]),$(ps[1])"
        @debug "max(xks): $(maximum(xks)), qs,ps[end]: $(qs[end]),$(ps[end])"

        ks = sind.(θs)
        @debug "ks before lookup: $ks"
        kis = [searchsortedlast(qs, k) for k in ks] |> collect
        ks = qs[kis]
        @debug "ks after lookup: $ks"
        θs = asind.(ks)

        # To access the the Fourier transform of the surface integral I(γ, q)
        # we must access the pattern ζ(p-q), so make a reverse index range for q
        # This works only if qs is symmetric (-Q/2, ..., 0, ..., Q/2)
        rev_qis = reverse(eachindex(qs))
        sFys_pqidxs = [i + qj - 1 for i in eachindex(ps), qj in rev_qis]
        @debug "sFys_pqidxs > Nx: \n$(findall(sFys_pqidxs .> Nx))"
        sFys_pqidxs[end, 1] = 1 # Alias with Nyquist frequency, avoid OOB error in sFys
        @debug "sFys_pqidxs > Nx: \n$(findall(sFys_pqidxs .> Nx))"
        if !all(qs[rev_qis] .≈ -qs)
            @warn "rev_qis reversed not matching qs: error=$(qs[rev_qis] .+ qs[rev_qis])"
        end

        rev_kis = rev_qis[kis]
        if !all(qs[rev_kis] .≈ -qs[kis])
            @warn "kis reversed not matching ks: error=$(qs[rev_kis] .+ qs[kis])"
        end

        seed = seed < 0 ? rand(0:typemax(Int)) : seed
        FFTplan = plan_fft!(similar(xs, ComplexF64))
        IFFTplan = plan_ifft!(similar(xs, ComplexF64))
        @debug "Parameters"
        @debug "kis: $kis"
        @debug "rev_kis: $rev_kis"
        @debug "rev_qis: $rev_qis"
        @debug "surf: $surf"
        @debug "above: $above"
        @debug "below: $below"
        @debug "seed: $seed"
        @debug "typeof(FFTplan): $(typeof(FFTplan))"
        @debug "typeof(IFFTplan): $(typeof(IFFTplan))"
        new{ST,AT,BT}(
            lambda, Lx, Nx, Nq, Ni, dx, dq,
            xs, xks, ps, qs,
            θs, ks, kis, sFys_pqidxs,
            surf, above, below,
            seed, Xoshiro(seed), FFTplan, IFFTplan,
        )
    end
end

function ParametersConfig(params::Parameters{ST,AT,BT}) where {ST<:RandomSurface,AT<:Material,BT<:Material}
    lambda = params.lambda
    k0 = 2π / lambda
    Lx = params.Lx / (k0 * lambda)
    Nx = params.Nx
    θs = params.θs
    Ni = params.Ni
    surf = scale(params.surf, 1.0 / k0)
    above = params.above
    below = params.below
    seed = params.seed
    return ParametersConfig(
        lambda=lambda,
        Lx=Lx,
        Nx=Nx,
        θs=θs,
        Ni=Ni,
        surf=surf,
        above=above,
        below=below,
        seed=seed
    )
end


function pretty(config::ParametersConfig)
    io = IOBuffer()
    println(io, "ParametersConfig(")
    for field in fieldnames(typeof(config))
        println(io, "\t$field=$(getfield(config, field)),")
    end
    println(io, ")")
    return String(take!(io))
end

function pretty(params::Parameters)
    config = ParametersConfig(params)
    return pretty(config)
end

Base.display(config::ParametersConfig) = println(pretty(config))
Base.display(params::Parameters) = println(pretty(params))

function save_parameters_config(file::String, params::ParametersConfig)
    file = endswith(file, ".conf") ? file : file * ".conf"
    open(file, "w") do io
        print(io, pretty(params))
    end
    return
end

function load_parameters_config(file::String)::ParametersConfig
    file = endswith(file, ".conf") ? file : file * ".conf"
    return eval(Meta.parse(read(file, String)))
end