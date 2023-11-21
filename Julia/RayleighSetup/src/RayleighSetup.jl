module RayleighSetup

#= Implements optical scattering under the reduced Rayleigh equations
# Effectively solves Maxwell's equations for a singularly polarized wave
# on a partially symmetric rough surface. Satisfying the boundary conditions

# The reduced Rayleigh equations are a set of coupled integral equations
# which assume that far field scattering conditions (singular direction, up/down in 1D)
# can be used all the way down to the rough surface boundary even though for strongly
# rough surfaces one can get multiple scattering events.
=#

using FFTW

export RayleighParams, SurfPreAlloc
export Polarization, s, p
export SurfType, flat, gaussian, singlebump
export show_params, gaussian_surface_gen, single_bump_gen
export gg, Wg
export test_rp, test_sp, test_rp_and_sp, c0

const c0 = 299_792_458.0 # [m/s], Speed of light in vacuum

@enum Polarization p s
struct RayleighParams

    #= 
    These are all parameters which doesn't
    during different realisations of the surface
    =#
    # Planned Fourier transform of surface points
    # cFFTWPlan{T,K,inplace,dims,flags}
    FT::FFTW.cFFTWPlan{ComplexF64,-1,true,1,Tuple{Int64}}

    xs::Vector{Float64}  # Surface points
    ps::Vector{Float64}  # Scattered wave numbers
    qs::Vector{Float64}  # Scattered wave numbers
    wq::Vector{Float64}  # Weights of the quadrature in q_n

    ν::Polarization # Polarization [p, s]
    ε::ComplexF64 # Permittivity of the scattering medium
    μ::ComplexF64 # Permeability of the scattering medium
    λ::Float64   # wavelength in [μm]
    ω::Float64   # Angular frequency

    # Sizings
    Q::Float64   # Truncated wave number, q = (-∞,∞) -> q_n = (-Q_mult / 2, Q_mult / 2)
    Δq::Float64  # Wave number spacing [1/m]
    Lx::Float64  # Surface length [m]
    Δx::Float64  # Surface point spacing [m]
    δ::Float64   # RMS height of surface profile function
    a::Float64   # Autocorrelation length

    Ni::Int      # Order of surface power expansion
    Nq::Int      # Number of quadratures in q and p, sizeof points give +1
    Nx::Int      # Number of surfaces sections in xs

    # Constructor
    RayleighParams(; ν::Polarization=p, ε=2.25 + 1e-4im, μ=1.0, λ=600e-9, Q_mult=4, Nq=127, L=10.0e-6, Ni=10, δ=30e-9, a=100e-9) = begin
        #=
        All lengths are being scaled to ω/c
        =#

        K = 2π / λ
        ω = c0 * K
        # All variables scaled such that ω/c0 = 1
        @assert imag(ε * μ) != 0.0 "ε times μ must have a small imaginary component to avoid singularities, but is $ε."

        # Assertions and warnings
        @assert L / λ > 10.0 "Surface length must be much larger than the wavelength, L ≫ λ, but is L:$L and λ:$λ."
        @assert Q_mult > 2 "Q_mult must be greater than 2, but is $Q_mult. ¤ is recommended."
        @assert Nq > 2 "Nq must be greater than 2, but is $Nq."

        Nx = 2 * Nq

        # Q = Q_mult * ω / c0 # Truncated wave number
        Q = Q_mult # Truncated wave number, divided by ω / c0
        Δq = Q / Nq

        ps = -Q/2:Δq:Q/2
        qs = Q/2:-Δq:-Q/2

        wq = ones(size(qs))
        wq[1] = wq[end] = 3.0 / 8.0
        wq[2] = wq[end-1] = 7.0 / 6.0
        wq[3] = wq[end-2] = 23.0 / 24.0

        Lx = L * ω / c0 # Surface length, scaled up by ω / c0, since reciprocal space is scaled down by ω / c0
        Δx = Lx / Nx
        xs = -Lx/2:Δx:Lx/2

        δ_new = δ * ω / c0 # RMS height of surface profile function scaled
        a_new = a * ω / c0 # Autocorrelation length scaled

        new(plan_fft!(similar(xs, ComplexF64)),
            xs, ps, qs, wq,
            ν, ε, μ, λ, ω,
            Q, Δq, Lx, Δx, δ_new, a_new,
            Ni, Nq, Nx)
    end
end


Wg(x, a) = exp(-x^2 / a^2)
gg(k, a) = √π * a * exp(-(a * k / 2.0)^2)

function gaussian_surface_gen(rp::RayleighParams)

    xs = rp.xs
    N = length(xs)
    ks = -rp.Q:rp.Δq:rp.Q

    Z = rp.δ .* randn(Float64, N) |> complex
    F = Wg.(xs, rp.a)

    return rp.FT \ ((rp.FT * Z) .* sqrt.(rp.FT * F))
end

function single_bump_gen(rp::RayleighParams)
    # In this case, δ is taken as the peak height
    # and a is taken as the standard deviation width of the gaussian bump
    return rp.δ * exp.((-0.5 / rp.a^2) * rp.xs .^ 2)
end

@enum SurfType flat gaussian singlebump
struct SurfPreAlloc
    #=
        SurfPreAlloc is a struct which contains all
        preallocated variables which are used during
        the calculation of the surface integral, i.e. the members are mutated
    =#

    # Preallocated steps in the calculations for a given surface
    Mpq::Matrix{ComplexF64}  # Matrix of the Mpq coefficients (A)
    Npk::Vector{ComplexF64}  # Vector of the Npk coefficients (b)
    R::Vector{ComplexF64}    # Reflection coefficient
    Fys::Vector{ComplexF64}  # Fourier transform of surface heights, prealloc
    sFys::Vector{ComplexF64} # Shifted Fourier transform of surface heights, prealloc
    ys::Vector{Float64}  # Surface heights

    function SurfPreAlloc(rp::RayleighParams, surf_t::SurfType)
        xs = rp.xs

        if surf_t == singlebump
            ys = single_bump_gen(rp)
        elseif surf_t == gaussian
            ys = gaussian_surface_gen(rp) .|> real
        elseif surf_t == flat
            ys = zeros(size(xs)) # Flat surface
        end

        Fys = similar(ys, ComplexF64)
        sFys = similar(Fys)
        Mpq = Matrix{ComplexF64}(undef, rp.Nq + 1, rp.Nq + 1)
        Npk = Vector{ComplexF64}(undef, rp.Nq + 1)
        R = similar(Npk)

        new(Mpq, Npk, R, Fys, sFys, ys)
    end
end

function show_params(rp::RayleighParams)
    names = fieldnames(typeof(rp))
    for name in names
        field = getfield(rp, name)
        Base.print("$name:\t")
        display(field)
    end
end

function show_params(sp::SurfPreAlloc)
    names = fieldnames(typeof(sp))
    for name in names
        field = getfield(sp, name)
        Base.print("$name:\t")
        display(field)
    end
end

function test_rp()
    Nq = 2^11 # = 2048
    rp = RayleighParams(;
        ν=p,
        Nq=Nq,
        ε=2.25,
        L=10.0e-6,
        Q_mult=4,
        Ni=10
    )
    display("RayleighParams test yielded:")
    show_params(rp)
    return rp
end

function test_sp(rp::RayleighParams, surf_t::SurfType=flat)
    sp = SurfPreAlloc(rp, surf_t)
    display("SurfPreAlloc test yielded:")
    show_params(sp)
    return sp
end

function test_rp_and_sp()
    rp = test_rp()
    sp_flat = test_sp(rp)
    sp_gauss = test_sp(rp, gauss)
    return rp, sp_flat, sp_gauss
end


end # module RayleighSetup