#= Implements optical scattering under the reduced Rayleigh equations
# Effectively solves Maxwell's equations for a singularly polarized wave
# on a partially symmetric rough surface. Satisfying the boundary conditions

# The reduced Rayleigh equations are a set of coupled integral equations
# which assume that far field scattering conditions (singular direction, up/down in 1D)
# can be used all the way down to the rough surface boundary even though for strongly
# rough surfaces one can get multiple scattering events.
=#

const c0 = 299_792_458.0 # [m/s], Speed of light in vacuum

@enum Polarization p s
@enum SurfType flat gaussian singlebump rect

polarization_from_string(str::String) = str == "p" ? p :
                                        str == "s" ? s :
                                        error("Unknown polarization: $str")

surftype_from_string(str::String) = str == "flat" ? flat :
                                    str == "gaussian" ? gaussian :
                                    str == "singlebump" ? singlebump :
                                    str == "rect" ? rect :
                                    error("Unknown surface type: $str")

struct Surface
    surf_t::SurfType
    surf_params::Vector{Float64}

    function Surface(surf_t::SurfType, surf_params::Vector{Float64})
        if surf_t == rect
            @assert length(surf_params) == 3 "rect requires 3 parameters, but got $(length(surf_params))"
        elseif surf_t == gaussian || surf_t == singlebump
            @assert length(surf_params) == 2 "$(surf_t) requires 2 parameters, but got $(length(surf_params))"
        else
            @assert length(surf_params) == 0 "$(surf_t) requires 0 parameters, but got $(length(surf_params))"
        end
        new(surf_t, surf_params)
    end
end

function scale(surf::Surface, scaling::Float64)::Surface
    # Scales the surface parameters in surf_params
    # by the scaling factor scaling

    surf_t = surf.surf_t

    if surf_t == rect
        δ = surf.surf_params[1] * scaling
        km = surf.surf_params[2]
        kp = surf.surf_params[3]
        return Surface(surf_t, [δ, km, kp])
    elseif surf_t == gaussian || surf_t == singlebump
        δ = surf.surf_params[1] * scaling
        a = surf.surf_params[2] * scaling
        return Surface(surf_t, [δ, a])
    else
        return Surface(surf_t, Vector{Float64}())
    end
end

as_string(surf::Surface)::String = (
    "surf_t=$(surf.surf_t)\n" *
    "surf_params=$(surf.surf_params)\n"
)

struct RayleighParams

    #= 
    These are all parameters which doesn't change
    during different realisations of the surface
    =#

    # Planned Fourier transform of surface points
    # cFFTWPlan{T,K,inplace,dims,flags}
    FT::FFTW.cFFTWPlan{ComplexF64,-1,true,1,Tuple{Int64}}

    xs::Vector{Float64}  # Surface points
    ps::Vector{Float64}  # Scattered wave numbers
    qs::Vector{Float64}  # Scattered wave numbers
    ks::Vector{Float64}  # Incident wave numbers
    kis::Vector{Int64}  # Indices of k in qs

    ν::Polarization # Polarization [p, s]
    ε::ComplexF64 # Permittivity of the scattering medium
    μ::ComplexF64 # Permeability of the scattering medium
    λ::Float64   # wavelength in [μm]
    ω::Float64   # Angular frequency

    # Sizings
    Q::Int64   # Truncated wave number, q = (-∞,∞) -> q_n = (-Q_mult / 2, Q_mult / 2)
    Δq::Float64  # Wave number spacing [1/m]
    Lx::Float64  # Surface length [m]
    Δx::Float64  # Surface point spacing [m]
    scaled::String # Scaling factor

    Ni::Int      # Order of surface power expansion
    Nq::Int      # Number of quadratures in q and p, sizeof points give +1
    Nx::Int      # Number of surfaces sections in xs

    surf::Surface # Parametrization of the surface

    # Constructor
    function RayleighParams(; ν::Polarization=p, ε=2.25, μ=1.0, λ=600e-9, Q=4, Nq=127, ks=[sind(10.0)], L=10.0e-6, Ni=10, surf::Surface, scaled="")::RayleighParams
        #=
        All lengths are being normalized to 'scaling'
        Assumes the surface lengths has already been scaled
        =#

        K = 2π / λ
        ω = c0 * K

        if imag(ε * μ) == 0.0
            println("εμ is purely real, adding a 1e-4i to ε to avoid singularities")
            ε += 1e-4im # Add a small imaginary part to avoid singularities
        end

        if scaled == ""
            surf = scale(surf, ω / c0)
            Lx = L * ω / c0
            scaled = "ω/c"
        else
            Lx = L
        end

        # Assertions and warnings
        @assert L / λ > 10.0 "Surface length must be much larger than the wavelength, L ≫ λ, but is L:$L and λ:$λ."
        @assert Q > 2 "Q must be greater than 2, but is $Q. ¤ is recommended."
        @assert Nq > 2 "Nq must be greater than 2, but is $Nq."


        Nx = 2 * Nq
        Δq = Q / Nq

        ps = -Q/2:Δq:Q/2
        qs = Q/2:-Δq:-Q/2

        Δx = Lx / Nx
        xs = -Lx/2:Δx:Lx/2

        kis = [searchsortedfirst(qs, k, rev=true) for k in ks]

        new(plan_fft!(similar(xs, ComplexF64)),
            xs, ps, qs, ks, kis,
            ν, ε, μ, λ, ω,
            Q, Δq, Lx, Δx, scaled,
            Ni, Nq, Nx, surf)
    end
end

function get_angles(rp::RayleighParams)::Vector{Float64}
    # Returns the angles of the incindent wave vectors
    # in the RayleighParams struct
    return asind.(rp.qs[rp.kis])
end

function as_string(rp::RayleighParams)::String
    # Returns a string representation of the RayleighParams struct
    return (
        "ν=$(rp.ν)\n" *
        "ε=$(rp.ε)\n" *
        "μ=$(rp.μ)\n" *
        "λ=$(rp.λ)\n" *
        "ω=$(rp.ω)\n" *
        "Q=$(rp.Q)\n" *
        "ks=$(rp.ks)\n" *
        "scaled=$(rp.scaled)\n" *
        "Nq=$(rp.Nq)\n" *
        "Nx=$(rp.Nx)\n" *
        "Lx=$(rp.Lx)\n" *
        "Ni=$(rp.Ni)\n" *
        as_string(rp.surf)
    )
end

function parse(::Type{RayleighParams}, str::String)::RayleighParams

    # Parse a string into a RayleighParams struct
    # The string must be formatted as the output of
    # string(RayleighParams)

    # Parse the string into a dictionary
    dict = Dict{String,String}()
    for line in split(str, '\n')
        if line != ""
            key, value = split(line, '=')
            dict[key] = value
        end
    end

    surf_params = dict["surf_params"]
    display(surf_params)
    if surf_params == "Float64[]"
        surf_params = Vector{Float64}()
    else
        surf_params = strip(surf_params, ['[', ']'])
        surf_params = parse.(Float64, split(surf_params, ','))
    end

    ks = dict["ks"]
    ks = strip(ks, ['[', ']'])
    # Parse the dictionary into a RayleighParams struct
    return RayleighParams(
        ν=polarization_from_string(dict["ν"]),
        ε=parse(ComplexF64, dict["ε"]),
        μ=parse(ComplexF64, dict["μ"]),
        λ=parse(Float64, dict["λ"]),
        Q=parse(Int64, dict["Q"]),
        ks=parse.(Float64, split(ks, ',')),
        Nq=parse(Int64, dict["Nq"]),
        L=parse(Float64, dict["Lx"]),
        Ni=parse(Int64, dict["Ni"]),
        scaled=dict["scaled"],
        surf=Surface(
            surftype_from_string(dict["surf_t"]),
            surf_params
        )
    )
end

function RayleighParams(filename::String)::RayleighParams
    # Loads a RayleighParams struct from an input file
    open(filename, "r") do io
        return parse(RayleighParams, read(io, String))
    end
end


#### O'Donnel rectangular correlation function for surface generation
H(x::Float64)::Float64 = (x > 0.0 ? 1.0 : 0.0) # Heaviside step function

Wr(x, km, kp) = (sin(kp * x) - sin(km * x)) / ((kp - km) * x) # Handle when x = 0 in rect_gen!

gr(k, km, kp) = π / (kp - km) * (H(kp - k) * H(k - km) + H(kp - k) * H(-k - km))

function rect_gen!(ys::Vector{Float64}, xs::Vector{Float64}, FT::FFTW.cFFTWPlan{ComplexF64,-1,true,1,Tuple{Int64}},
    δ::Float64, km::Float64, kp::Float64)::Nothing
    # Overwrites ys with a surface defined by the
    # rectangular West-O'Donnel correlation function
    # δ is the RMS height in unit of c/ω
    Z = δ .* randn(Float64, length(ys)) |> complex

    # xs → 0.0 is a special case, so we handle it separately
    zr = ceil(Int64, length(xs) / 2)

    corr = Vector{ComplexF64}(undef, length(xs))
    corr[zr] = 1.0 # From Taylor series around x = 0
    corr[zr+1:end] .= Wr.(xs[zr+1:end], km, kp)
    corr[1:zr-1] .= Wr.(xs[1:zr-1], km, kp)

    ys .= FT \ ((FT * Z) .* sqrt.(FT * corr)) .|> real
    return nothing
end

#### Gaussian correlation function for surface generation
Wg(x, a) = exp(-x^2 / a^2)
gg(k, a) = √π * a * exp(-(a * k / 2.0)^2)

function gaussian_gen!(ys::Vector{Float64}, xs::Vector{Float64}, FT::FFTW.cFFTWPlan{ComplexF64,-1,true,1,Tuple{Int64}},
    δ::Float64, a::Float64)::Nothing
    # Overwrites ys with a surface defined by
    # the Gaussian correlation function

    Z = δ .* randn(Float64, length(ys)) |> complex
    ys .= FT \ ((FT * Z) .* sqrt.(FT * Wg.(xs, a))) .|> real
    return nothing
end

#### Single bump correlation function for surface generation
function single_bump_gen!(ys::Vector{Float64}, xs::Vector{Float64},
    δ::Float64, a::Float64)::Nothing
    # In this case, δ is taken as the peak height
    # and a is taken as the RMS width of the bump
    ys .= δ * exp.((-0.5 / a^2) * xs .^ 2)
    return nothing
end

#### Flat correlation function for surface generation
function flat_gen!(ys::Vector{Float64}, xs::Vector{Float64})::Nothing
    ys .= zeros(length(xs))
    return nothing
end

struct SurfPreAlloc
    #=
        SurfPreAlloc is a struct which contains all
        preallocated variables which are used during
        the calculation of the surface integral, i.e. the members are mutated
        Note that all members are uninitialized, and must be initialized after making the struct
    =#

    # Preallocated steps in the calculations for a given surface
    Mpq::Matrix{ComplexF64}  # Matrix of the Mpq coefficients (A)
    Npk::Matrix{ComplexF64}  # Vectors (for each k) of the Npk coefficients (b)
    Fys::Vector{ComplexF64}  # Fourier transform of surface heights, prealloc
    sFys::Vector{ComplexF64} # Shifted Fourier transform of surface heights, prealloc
    ys::Vector{Float64}  # Surface heights

    function SurfPreAlloc(rp::RayleighParams)
        xs = rp.xs

        ys = similar(xs, Float64)
        Fys = similar(ys, ComplexF64)
        sFys = similar(Fys)
        Mpq = Matrix{ComplexF64}(undef, rp.Nq + 1, rp.Nq + 1)
        Npk = Matrix{ComplexF64}(undef, rp.Nq + 1, length(rp.kis))

        new(Mpq, Npk, Fys, sFys, ys)
    end
end

function surfgen_func_selector!(rp::RayleighParams)::Function
    # Returns a function which overwrites sp.ys for each call
    surf_t = rp.surf.surf_t
    if surf_t == rect
        # Qs = -rp.Q:rp.Δq:rp.Q |> collect
        return (ys) -> rect_gen!(ys, rp.xs, rp.FT, rp.surf.surf_params...)
    elseif surf_t == gaussian
        # Qs = -rp.Q:rp.Δq:rp.Q |> collect
        return (ys) -> gaussian_gen!(ys, rp.xs, rp.FT, rp.surf.surf_params...)
    elseif surf_t == singlebump
        return (ys) -> single_bump_gen!(ys, rp.xs, rp.surf.surf_params...)
    else
        return (ys) -> flat_gen!(ys, rp.xs)
    end
end

#### Save, load and make solver config
function save_solver_config(file::String, rp::RayleighParams)::Nothing
    # Saves the RayleighParams struct to a file
    # in the data folder
    open(file, "w") do io
        write(io, as_string(rp))
    end
    return nothing
end

function load_solver_config(file::String)::Tuple{RayleighParams,SurfPreAlloc,Function}
    # Loads a complete solver setup from a config file
    # Config file should ideally be generated with save_solver_config
    # Returns a RayleighParams struct, a SurfPreAlloc struct and a surface generator function
    rp = RayleighParams(file)
    sp = SurfPreAlloc(rp)
    surfgen_func! = surfgen_func_selector!(rp)

    return rp, sp, surfgen_func!
end

function make_solver_config()::Tuple{RayleighParams,SurfPreAlloc,Function}
    print("Input for solver parameters in RayleighParams struct\n")
    print("Polarization [p|s] (=p): ")
    in = readline()
    ν = polarization_from_string(in == "" ? "p" : in)

    print("ε [X+Yim] (=2.25): ")
    in = readline()
    ε = parse(ComplexF64, in == "" ? "2.25" : in)

    print("μ [X+Yim] (=1.0): ")
    in = readline()
    μ = parse(ComplexF64, in == "" ? "1.0" : in)

    print("λ [nm] (=632.8): ")
    in = readline()
    λ = parse(Float64, in == "" ? "600" : in) * 1e-9

    print("Q [multiple of ω/c] (=4): ")
    in = readline()
    Q = parse(Int64, in == "" ? "4" : in)

    print("θs [list of deg OR \"fresnel\"] (=0,10,20): ")
    in = readline()
    if in == "fresnel"
        θs = 0.0:0.5:90.0
    else
        θs = parse.(Float64, split(in == "" ? "0,10,20" : in, ','))
    end

    print("Nq (=1024): ")
    in = readline()
    Nq = parse(Int64, in == "" ? "1024" : in)

    print("L [multiple of λ] (=100): ")
    in = readline()
    L = parse(Float64, in == "" ? "100" : in) * λ

    print("Ni (=10): ")
    in = readline()
    Ni = parse(Int64, in == "" ? "10" : in)

    print("Surface type [flat|gaussian|singlebump|rect] (=gaussian): ")
    in = readline()
    surf_t = surftype_from_string(in == "" ? "gaussian" : in)

    if surf_t == rect
        print("δ [nm] (=5.0): ")
        in = readline()
        δ = parse(Float64, in == "" ? "5.0" : in) * λ

        print("km [scaled to ω/c] (=0.8): ")
        in = readline()
        km = parse(Float64, in == "" ? "0.8" : in)

        print("kp [scaled to ω/c] (=1.2): ")
        in = readline()
        kp = parse(Float64, in == "" ? "1.2" : in)

        surf_params = [δ, km, kp]
    elseif surf_t == gaussian || surf_t == singlebump
        print("δ [nm] (=5.0): ")
        in = readline()
        δ = parse(Float64, in == "" ? "5.0" : in) * λ

        print("a [multiple of λ] (=1.0): ")
        in = readline()
        a = parse(Float64, in == "" ? "1.0" : in) * λ

        surf_params = [δ, a]
    else
        surf_params = Vector{Float64}()
    end

    rp = RayleighParams(
        ν=ν,
        ε=ε,
        μ=μ,
        λ=λ,
        Q=Q,
        ks=sind.(θs),
        Nq=Nq,
        L=L,
        Ni=Ni,
        surf=Surface(surf_t, surf_params)
    )
    sp = SurfPreAlloc(rp)
    surfgen_func! = surfgen_func_selector!(rp)

    print("Save config to file? [y|n] (=n): ")
    in = readline()
    if in == "y"
        print("Filename: ")
        in = readline()
        save_solver_config(joinpath("input", in), rp)
    end

    return rp, sp, surfgen_func!
end