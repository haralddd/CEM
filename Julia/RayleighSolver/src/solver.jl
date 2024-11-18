
struct SimOutput
    R::Matrix{ComplexF64}
    R²::Matrix{Float64}
    σ²::Matrix{Float64}
    κ::Matrix{Float64}

    function SimOutput(spa::SimParams)
        # Reduced q-vector indices
        R = zeros(ComplexF64, (length(spa.qs), length(spa.ks)))
        R² = zeros(Float64, (length(spa.qs), length(spa.ks)))
        σ² = zeros(Float64, (length(spa.qs), length(spa.ks)))
        κ = zeros(Float64, (length(spa.qs), length(spa.ks)))
        return new(R, R², σ², κ)
    end
end

struct SolverData{SimParamsType}
    spa::SimParamsType
    sp::SimPrealloc
    out_p::SimOutput
    out_s::SimOutput
    iters::Int64

    function SolverData(spa::SimParamsType, iters::Int64=1000) where {SimParamsType}
        out_p = SimOutput(spa)
        @debug "SolverData init: SimOutput P-polarization complete"
        out_s = SimOutput(spa)
        @debug "SolverData init: SimOutput S-polarization complete"

        sp = SimPrealloc(spa)
        @debug "SolverData init: SimPrealloc complete"

        return new{SimParamsType}(spa, sp, out_p, out_s, iters)
    end
end

function precompute!(data::SolverData{SPA})::Nothing where {SPA}
    precompute!(data.sp, data.spa)
    validate(data.sp)
    return nothing
end

"Iteratively update observable like so: ⟨A⟩ₙ = (n-1)/n ⟨A⟩ₙ₋₁ + Aₙ/n"
function observe(obs, x, n)
    return obs + (x-obs)/n
end

"Updates all observables in out::SimOutput"
function observe!(out::SimOutput, A, n)
    R = out.R
    R² = out.R²
    σ² = out.σ²
    κ = out.κ

    @inbounds for I in eachindex(R)
        R[I] = observe(R[I], A[I], n)
    end

    @inbounds for I in eachindex(R²)
        R²[I] = observe(R²[I], abs2(A[I]), n)
    end

    @inbounds for I in eachindex(σ²)
        var = abs2(A[I] - R[I])
        σ²[I] = observe(σ²[I], var, n)
        κ[I] = observe(κ[I], var^2, n)
    end
end

function observe!(data::SolverData, n::Int)
    R_p = data.sp.p_data.Npk
    R_s = data.sp.s_data.Npk
    observe!(data.out_p, R_p, n)
    observe!(data.out_s, R_s, n)
end

"""
    function solve_single!(data::SolverData)::Nothing

Calculates the preallocated surface integral.
Matrix Mpqn is preallocated and contains the invariant parts of the Mpq matrix.
Matrix Npkn is preallocated and contains the invariant parts of the Npk vector.
Stores the resulting reflection coefficients in sp.Npk
Overwrites sp.Mpq with the LU factorization in the linear solution process.

# Arguments:
- `data`: [`SolverData`](@ref) - Contains the parameters, preallocated steps, output
"""
function solve_single!(data::SolverData)::Nothing

    spa = data.spa
    sp = data.sp

    FFT = spa.FFT

    sFys_pqidxs = spa.sFys_pqidxs
    kis = spa.kis

    ys = sp.ys
    Fys = sp.Fys
    sFys = sp.sFys

    pd = sp.p_data
    sd = sp.s_data

    pd.Mpq .= 0.0
    pd.Npk .= 0.0

    sd.Mpq .= 0.0
    sd.Npk .= 0.0
    
    @inbounds for n in reverse(axes(pd.Mpqn, 3)) # Reverse because prefactors vanish at higher powers of ´n´
        for i in eachindex(Fys)
            Fys[i] = ys[i]^(n - 1)
        end

        FFT * Fys # In place FFT
        fftshift!(sFys, Fys)

        for j in axes(pd.Mpq, 2)
            for i in axes(pd.Mpq, 1)
                idx = sFys_pqidxs[i, j]
                pd.Mpq[i, j] += pd.Mpqn[i, j, n] * sFys[idx]
                sd.Mpq[i, j] += sd.Mpqn[i, j, n] * sFys[idx]
            end
        end
        for (j, kj) in enumerate(kis)
            for i in axes(pd.Npk, 1)
                idx = sFys_pqidxs[i, kj]
                pd.Npk[i, j] += pd.Npkn[i, j, n] * sFys[idx]
                sd.Npk[i, j] += sd.Npkn[i, j, n] * sFys[idx]
            end
        end
    end

    A = lu!(pd.Mpq)
    @inbounds for j in axes(pd.Npk, 2)
        b = @view pd.Npk[:, j]
        ldiv!(A, b)
        # ldiv!(A, pd.Npk[:, j])
    end

    A = lu!(sd.Mpq)
    @inbounds for j in axes(sd.Npk, 2)
        b = @view sd.Npk[:, j]
        ldiv!(A, b)
        # ldiv!(A, pd.Npk[:, j])
    end

    return nothing
end


function solve_MDRC!(data::SolverData)
    # Solve the Rayleigh problem for a given SimParams struct
    # sp is the preallocated surface struct, containing the preallocated output, this is mutated in place
    # spa is the SimParams struct, containing the constant parameters for the calculation
    # surf_generator! is a function that generates the surface sp.ys for each iteration
    # N_ens is the number of ensemble averages to perform
    # returns the coherent and incoherent MDRC

    @info "precompute!:"
    @time precompute!(data)

    thr_sz = Threads.nthreads()
    LinearAlgebra.BLAS.set_num_threads(thr_sz)

    sp = data.sp
    spa = data.spa
    iters = data.iters

    @info "mainloop:"

    for n in ProgressBar(1:iters)
    # @time for n in 1:iters
        generate_surface!(sp, spa)
        solve_single!(data)
        observe!(data, n)
    end

    return nothing
end

# function get_coh_inc(out::SimOutput, θss::Vector{Float64}, θis::Vector{Float64}, Lx::Float64, k0::Float64)
#     coh = Matrix{Float64}(undef, (length(θss), length(θis)))
#     inc = Matrix{Float64}(undef, (length(θss), length(θis)))
#     for (j, θi) in enumerate(θis)
#         for (i, θs) in enumerate(θss)
#             C = alp^2 / (2π*Lx*cos(θi))
#             coh[i, j] = C * abs2(out.R[i, j])
#             inc[i, j] = C * out.R²[i, j] - coh[i, j]
#         end
#     end
#     return coh, inc
# end

function get_coh_inc(R, R2, qs::Vector{Float64}, ks::Vector{Float64})
    coh = Matrix{Float64}(undef, (length(qs), length(ks)))
    inc = Matrix{Float64}(undef, (length(qs), length(ks)))
    for (j, k) in enumerate(ks)
        for (i, q) in enumerate(qs)
            C = abs2(alpha0(q)) / abs(alpha0(k))
            coh[i, j] = C * abs2(R[i, j])
            inc[i, j] = C * R2[i, j] - coh[i, j]
        end
    end
    return coh, inc
end

struct DataMDRC
    coh_p::Matrix{Float64}
    inc_p::Matrix{Float64}
    coh_s::Matrix{Float64}
    inc_s::Matrix{Float64}
    θis::Vector{Float64}
    θss::Vector{Float64}
end

"Returns coherent and incoherent MDRC calculated from ⟨R⟩ and ⟨R²⟩"
function calc_mdrc(data::SolverData)
    spa = data.spa
    qs = spa.qs
    ks = spa.ks

    mask = qs .>= -1.0 .&& qs .<= 1.0
    θss = asind.(qs[mask])
    θis = data.spa.θs

    @assert all(ks .≈ sind.(θis))
    Rp = data.out_p.R
    R2p = data.out_p.R²

    Rs = data.out_s.R
    R2s = data.out_s.R²

    coh_p, inc_p = get_coh_inc(Rp[mask, :], R2p[mask,:], qs[mask], ks)
    coh_s, inc_s = get_coh_inc(Rs[mask, :], R2s[mask, :], qs[mask], ks)

    return DataMDRC(coh_p, inc_p, coh_s, inc_s, θis, θss)
end
