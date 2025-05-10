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
    coh_p::Matrix{Float64}
    inc_p::Matrix{Float64}
    coh_s::Matrix{Float64}
    inc_s::Matrix{Float64}
    θtps::Vector{Float64}
    θtss::Vector{Float64}
    θtes::Vector{Float64}
    θtos::Vector{Float64}
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
            C = Lx/2π * real(alpha0(q)) * real(alpha0(q) / alpha0(k))
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
    below = params.below
    if below isa Uniaxial
        κpa = ν == :p ? below.eps_para : below.mu_para
        alpha_func = ν == :p ? (q -> alpha_p(q, below)) : (q -> alpha_s(q, below))
    elseif below isa Isotropic
        κpa = below.mu
        alpha_func = (q -> alpha(q, below))
    end
    for (j, k) in enumerate(ks)
        for (i, q) in enumerate(qs)
            a = alpha_func(q)
            C = Lx/2π * abs(real(a) * real(a / (κpa * alpha0(k))))
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
    θ0s = params.θs

    @assert all(ks .≈ sind.(θ0s))
    Rp = data.Rp.A[mask, :]
    R2p = data.Rp.A²[mask, :]

    Rs = data.Rs.A[mask, :]
    R2s = data.Rs.A²[mask, :]

    coh_p, inc_p = get_mdrc_coh_inc(Rp, R2p, qs[mask], ks, params)
    coh_s, inc_s = get_mdrc_coh_inc(Rs, R2s, qs[mask], ks, params)

    return MdrcPlotData(coh_p, inc_p, coh_s, inc_s, θss, θ0s)
end

"""
    θt(q::Float64, alpha::Float64)

Calculates the transmission angle in degrees from the wavenumber and alpha value.

# Arguments:
- `q`: Wavenumber component
- `alpha`: Alpha value for the material
"""
function θt(q::Float64, alpha::Float64)
    return atand(q, abs(alpha))
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
function calc_mdtc(data::SolverData)
    params = data.params
    qs = params.qs
    ks = params.ks
    below = data.params.below
    if below isa Uniaxial
        mu = below.mu_perp
        μεpe = mu * below.eps_perp
        μεpa = mu * below.eps_para
        ap_func = (q -> alpha_p(q, below))
        as_func = (q -> alpha_s(q, below))
    elseif below isa Isotropic
        mu = below.mu
        μεpe = mu * below.eps
        μεpa = mu * below.eps
        ap_func = (q -> alpha(q, below))
        as_func = ap_func
    end
    
    mask_p = qs .>= -real(sqrt(μεpe)) .&& qs .<= real(sqrt(μεpe))
    mask_s = qs .>= -real(sqrt(μεpa)) .&& qs .<= real(sqrt(μεpa))
    
    ap = ap_func.(qs[mask_p])
    as = as_func.(qs[mask_s])

    θtp = θt.(qs[mask_p], real.(ap))
    θts = θt.(qs[mask_s], real.(as))
    
    # Flip the angles if the real part is negative
    θtp[real(ap) .< 0.0] .= .-θtp[real(ap) .< 0.0]
    θts[real(as) .< 0.0] .= .-θts[real(as) .< 0.0]
    
    apk = ap_func.(ks)
    ask = as_func.(ks)
    
    θte = θt.(ks, real(apk))
    θto = θt.(ks, real(ask))

    # Flip the angles if the real part is negative
    θte[real(apk) .< 0.0] .= .-θte[real(apk) .< 0.0]
    θto[real(ask) .< 0.0] .= .-θto[real(ask) .< 0.0]

    Tp = data.Tp.A[mask_p, :]
    T2p = data.Tp.A²[mask_p, :]

    Ts = data.Ts.A[mask_s, :]
    T2s = data.Ts.A²[mask_s, :]

    coh_p, inc_p = get_mdtc_coh_inc(Tp, T2p, qs[mask_p], ks, params, :p)
    coh_s, inc_s = get_mdtc_coh_inc(Ts, T2s, qs[mask_s], ks, params, :s)

    return MdtcPlotData(coh_p, inc_p, coh_s, inc_s, θtp, θts, θte, θto)
end