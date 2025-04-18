"""
    assert_branch_cuts(qs, material)

Validates that the branch cuts are correctly defined for all wavenumbers in `qs`.
Ensures that the imaginary part of alpha is non-negative and the real part is non-negative.

# Arguments:
- `qs`: Vector of wavenumbers to check
- `material`: Material object containing the properties to calculate alpha
"""
function assert_branch_cuts(qs, material)
    for q in qs
        @assert imag(alpha(q, material)) >= 0.0
        @assert real(alpha(q, material)) >= 0.0
    end
end


"""
    solve_MDRC!(data::SolverData{_P}) where {_P}

Solves the Rayleigh problem for multiple surface realizations and computes ensemble averages.
Handles both multi-threaded and single-threaded computation depending on the system configuration.

# Arguments:
- `data`: [`SolverData`](@ref) - Contains the parameters, preallocated steps, and output structures
"""
function solve_MDRC!(data::SolverData{_P}) where {_P}


    # @info "Running tests:"
    # @info data.params.above
    # @info data.params.below
    # assert_branch_cuts(data.params.qs, data.params.above)
    # assert_branch_cuts(data.params.qs, data.params.below)
    # @info "Branch cuts OK"

    @info "Precomputing matrix elements:"
    precomputed = Precomputed(data.params)
    @time precompute!(precomputed, data.params)
    validate(precomputed)
    do_debug = get(ENV, "JULIA_DEBUG", "") != ""
    show_iter = get(ENV, "JULIA_SHOWITERS", "") == "true"
    @info "debug: $do_debug"
    @info "show iters: $show_iter"

    Nthreads = Threads.nthreads()
    LinearAlgebra.BLAS.set_num_threads(1)
    Niters = data.iters

    enable_threads = Nthreads > 1
    iter = 0

    if Niters > Nthreads && enable_threads
        @info "Mainloop: Running $Niters iters on $Nthreads threads..."
        
        allocs = [Preallocated(data.params) for _ in 1:Nthreads]
        lk = ReentrantLock()

        @time Threads.@threads :static for _ in (show_iter ? ProgressBar(1:Niters) : 1:Niters)
            tidx = Threads.threadid()
            _alloc = allocs[tidx]

            generate_surface!(_alloc, data.params)
            solve_single!(_alloc, precomputed, data)
            @lock lk begin
                iter += 1
                observe!(data.P_res, _alloc.PNpk, iter)
                observe!(data.S_res, _alloc.SNpk, iter)
                if do_debug
                    P_energy, S_energy = energy_conservation(abs2(_alloc.PNpk), abs2(_alloc.SNpk), data.params)
                    @debug "P_energy: $P_energy"
                    @debug "S_energy: $S_energy"
                end
            end # lock
        end
    else
        @info "Mainloop: Running $Niters iters sequentially."
        alloc = Preallocated(data.params)
        @time for n in (show_iter ? ProgressBar(1:Niters) : 1:Niters)
            generate_surface!(alloc, data.params)
            solve_single!(alloc, precomputed, data)
            observe!(data.P_res, alloc.PNpk, n)
            observe!(data.S_res, alloc.SNpk, n)
            iter += 1
            if do_debug
                P_energy, S_energy = energy_conservation(abs2(alloc.PNpk), abs2(alloc.SNpk), data.params)
                @debug "P_energy: $P_energy"
                @debug "S_energy: $S_energy"
            end
        end
    end


    return nothing
end

"""
    MdrcPlotData

Structure containing coherent and incoherent mean differential reflection coefficients for plotting.

# Fields:
- `coh_p`: Matrix of coherent P-polarization MDRC values
- `inc_p`: Matrix of incoherent P-polarization MDRC values
- `coh_s`: Matrix of coherent S-polarization MDRC values
- `inc_s`: Matrix of incoherent S-polarization MDRC values
- `θs`: Vector of scattering angles (degrees)
- `θ0s`: Vector of incident angles (degrees)
"""
struct MdrcPlotData
    coh_p::Matrix{Float64}
    inc_p::Matrix{Float64}
    coh_s::Matrix{Float64}
    inc_s::Matrix{Float64}
    θs::Vector{Float64}
    θ0s::Vector{Float64}
end

"""
    MdtcPlotData

Structure containing coherent and incoherent mean differential transmission coefficients for plotting.

# Fields:
- `coh`: Matrix of coherent MDTC values
- `inc`: Matrix of incoherent MDTC values
- `θs`: Vector of scattering angles (degrees)
- `θtes`: Vector of extraordinary transmitted angles (degrees)
- `θtos`: Vector of ordinary transmitted angles (degrees)
"""
struct MdtcPlotData
    coh::Matrix{Float64}
    inc::Matrix{Float64}
    θs::Vector{Float64}
    θtes::Vector{Float64}
    θtos::Vector{Float64}
end

"""
    θt(q::Float64, alpha::Float64)

Calculates the transmission angle in degrees from the wavenumber and alpha value.

# Arguments:
- `q`: Wavenumber component
- `alpha`: Alpha value for the material
"""
function θt(q::Float64, alpha::Float64)
    return atand(q, alpha)
end

"""
    get_mdrc_coh_inc(R, R2, qs::Vector{Float64}, ks::Vector{Float64}, params::Parameters)

Calculates the coherent and incoherent mean differential reflection coefficients.

# Arguments:
- `R`: Matrix of ensemble-averaged reflection amplitudes
- `R2`: Matrix of ensemble-averaged squared reflection amplitudes
- `qs`: Vector of scattered wavenumbers
- `ks`: Vector of incident wavenumbers
- `params`: [`Parameters`](@ref) containing simulation parameters

# Returns:
- Tuple of (coherent MDRC, incoherent MDRC) matrices
"""
function get_mdrc_coh_inc(R, R2, qs::Vector{Float64}, ks::Vector{Float64}, params::Parameters)
    Lx = params.Lx
    λ = params.lambda

    coh = Matrix{Float64}(undef, (length(qs), length(ks)))
    inc = Matrix{Float64}(undef, (length(qs), length(ks)))
    for (j, k) in enumerate(ks)
        for (i, q) in enumerate(qs)
            C = Lx/2π * alpha0(q)^2 / alpha0(k)
            coh[i, j] = C * abs2(R[i, j])
            inc[i, j] = C * R2[i, j] - coh[i, j]
        end
    end
    return coh, inc
end

"""
    get_mdtc_coh_inc(T, T2, qs::Vector{Float64}, ks::Vector{Float64}, params::Parameters, ν::Symbol=:p)

Calculates the coherent and incoherent mean differential transmission coefficients.

# Arguments:
- `T`: Matrix of ensemble-averaged transmission amplitudes
- `T2`: Matrix of ensemble-averaged squared transmission amplitudes
- `qs`: Vector of scattered wavenumbers
- `ks`: Vector of incident wavenumbers
- `params`: [`Parameters`](@ref) containing simulation parameters
- `ν`: Polarization, either `:p` or `:s`

# Returns:
- Tuple of (coherent MDTC, incoherent MDTC) matrices
"""
function get_mdtc_coh_inc(T, T2, qs::Vector{Float64}, ks::Vector{Float64}, params::Parameters, ν::Symbol=:p)
    Lx = params.Lx
    λ = params.lambda

    coh = Matrix{Float64}(undef, (length(qs), length(ks)))
    inc = Matrix{Float64}(undef, (length(qs), length(ks)))
    κpa = ν == :p ? params.below.eps_para : params.below.mu_para
    κpe = ν == :p ? params.below.eps_perp : params.below.mu_perp
    A_val = A(params.below)
    μεpe = params.below.mu_perp * params.below.eps_perp
    μεpa = params.below.mu_para * params.below.eps_para
    for (j, k) in enumerate(ks)
        for (i, q) in enumerate(qs)
            a = ν == :p ? alpha_p(q, A_val, μεpa) : alpha_s(q, μεpa)
            C = Lx/2π * real(a^2 / (κpa * alpha0(k)))
            coh[i, j] = C * abs2(T[i, j])
            inc[i, j] = C * T2[i, j] - coh[i, j]
        end
    end
    return coh, inc
end


"""
    calc_mdrc(data::SolverData)

Returns coherent and incoherent MDRC calculated from ⟨A⟩ and ⟨A²⟩.

# Arguments:
- `data`: [`SolverData`](@ref) containing the simulation results

# Returns:
- [`MdrcPlotData`](@ref) structure with coherent and incoherent reflection coefficients
"""
function calc_mdrc(data::SolverData)
    params = data.params
    qs = params.qs
    ks = params.ks

    mask = qs .>= -1.0 .&& qs .<= 1.0
    θss = asind.(qs[mask])
    θ0s = data.params.θs

    @assert all(ks .≈ sind.(θ0s))
    Rp = data.P_res.A
    R2p = data.P_res.A²

    Rs = data.S_res.A
    R2s = data.S_res.A²

    coh_p, inc_p = get_mdrc_coh_inc(Rp[mask, :], R2p[mask, :], qs[mask], ks, params)
    coh_s, inc_s = get_mdrc_coh_inc(Rs[mask, :], R2s[mask, :], qs[mask], ks, params)

    return MdrcPlotData(coh_p, inc_p, coh_s, inc_s, θss, θ0s)
end

"""
    calc_mdrc(data::SolverData{Parameters{_S,Vacuum,Uniaxial}}) where _S

Specialized implementation of calc_mdrc for Vacuum-Uniaxial interface.
Returns coherent and incoherent MDRC calculated from ⟨A⟩ and ⟨A²⟩.

# Arguments:
- `data`: [`SolverData`](@ref) containing the simulation results for a Vacuum-Uniaxial interface

# Returns:
- [`MdrcPlotData`](@ref) structure with coherent and incoherent reflection coefficients
"""
function calc_mdrc(data::SolverData{Parameters{_S,Vacuum,Uniaxial}}) where _S
    params = data.params
    qs = params.qs
    ks = params.ks

    mask = qs .>= -1.0 .&& qs .<= 1.0
    θss = asind.(qs[mask])
    θ0s = data.params.θs

    @assert all(ks .≈ sind.(θ0s))
    Rp = get_A(data.P_res)
    R2p = get_A²(data.P_res)

    Rs = get_A(data.S_res)
    R2s = get_A²(data.S_res)

    coh_p, inc_p = get_mdrc_coh_inc(Rp[mask, :], R2p[mask, :], qs[mask], ks, params)
    coh_s, inc_s = get_mdrc_coh_inc(Rs[mask, :], R2s[mask, :], qs[mask], ks, params)

    return MdrcPlotData(coh_p, inc_p, coh_s, inc_s, θss, θ0s)
end

"""
    calc_mdtc(data::SolverData{Parameters{_S,Vacuum,Uniaxial}}) where _S

Specialized implementation of calc_mdtc for Vacuum-Uniaxial interface.
Returns coherent and incoherent MDTC calculated from ⟨T⟩ and ⟨T²⟩.

# Arguments:
- `data`: [`SolverData`](@ref) containing the simulation results for a Vacuum-Uniaxial interface

# Returns:
- [`MdtcPlotData`](@ref) structure with coherent and incoherent transmission coefficients
"""
function calc_mdtc(data::SolverData{Parameters{_S,Vacuum,Uniaxial}}) where _S
    params = data.params
    qs = params.qs
    ks = params.ks
    below = data.params.below
    μεpe = below.mu_perp * below.eps_perp
    μεpa = below.mu_para * below.eps_para
    A = sqrt(μεpa/μεpe)

    # Transmission index ellipsoid is larger than reflection
    mask_p = qs .>= -real(sqrt(μεpe)) .&& qs .<= real(sqrt(μεpe))
    mask_s = qs .>= -real(sqrt(μεpa)) .&& qs .<= real(sqrt(μεpa))
    
    ap = real.(alpha_p.(qs[mask_p], A, μεpa))
    as = real.(alpha_s.(qs[mask_s], μεpa))
    
    θtp = θt.(qs[mask_p], ap)
    θts = θt.(qs[mask_s], as)
    
    θtes = [θt.(k, real(alpha_p(k, A, μεpa))) for k in ks]
    θtos = [θt.(k, real(alpha_s(k, μεpa))) for k in ks]

    Tp = get_T(data.P_res)
    T2p = get_T²(data.P_res)

    Ts = get_T(data.S_res)
    T2s = get_T²(data.S_res)

    coh_p, inc_p = get_mdtc_coh_inc(Tp[mask_p, :], T2p[mask_p, :], qs[mask_p], ks, params, :p)
    coh_s, inc_s = get_mdtc_coh_inc(Ts[mask_s, :], T2s[mask_s, :], qs[mask_s], ks, params, :s)

    return MdtcPlotData(coh_p, inc_p, θtp, θtes, θtos), MdtcPlotData(coh_s, inc_s, θts, θtes, θtos)
end