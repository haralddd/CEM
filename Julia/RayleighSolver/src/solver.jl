
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

    function SolverData(spa::SimParamsType, iters::Int64) where {SimParamsType}
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
function observe(observable, value, N)
    return observable*(N-1)/N + value/N
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

    rev_kis = spa.rev_kis
    rev_qis = spa.rev_qis

    ys = sp.ys
    Fys = sp.Fys
    sFys = sp.sFys

    pd = sp.p_data
    sd = sp.s_data

    pd.Mpq .= zero(ComplexF64)
    pd.Npk .= zero(ComplexF64)

    sd.Mpq .= zero(ComplexF64)
    sd.Npk .= zero(ComplexF64)

    @inbounds for n in reverse(axes(pd.Mpqn, 3)) # Reverse because prefactors vanish at higher powers of ´n´
        @inbounds for i in eachindex(Fys)
            Fys[i] = ys[i]^(n - 1)
        end

        FFT * Fys # In place FFT

        fftshift!(sFys, Fys)

        @inbounds for (j, qj) in enumerate(rev_qis)
            @inbounds for i in axes(pd.Mpq, 1)
                pd.Mpq[i, j] += pd.Mpqn[i, j, n] * sFys[i+qj-1]
                sd.Mpq[i, j] += sd.Mpqn[i, j, n] * sFys[i+qj-1]
            end
        end
        @inbounds for (j, kj) in enumerate(rev_kis)
            @inbounds for i in axes(pd.Npk, 1)
                pd.Npk[i, j] -= pd.Npkn[i, j, n] * sFys[i+kj-1]
                sd.Npk[i, j] -= sd.Npkn[i, j, n] * sFys[i+kj-1]
            end
        end
    end
    A_p = factorize(pd.Mpq)
    A_s = factorize(sd.Mpq)

    for i in eachindex(rev_kis)
        ldiv!(pd.Npk[:, i], pd.Mpq[:,:,i], pd.Npk)
        ldiv!(sd.Npk[:, i], sd.Mpq[:, :, i], sd.Npk)
    end
    
    # pd.Npk = pd.Mpq / pd.Npk
    # sd.Npk = sd.Mpq / pd.Npk


    # LinearAlgebra.LAPACK.gesv!(pd.Mpq, pd.Npk)
    # LinearAlgebra.LAPACK.gesv!(sd.Mpq, sd.Npk)

    return nothing
end


function solve_MDRC!(data::SolverData)
    # Solve the Rayleigh problem for a given SimParams struct
    # sp is the preallocated surface struct, containing the preallocated output, this is mutated in place
    # spa is the SimParams struct, containing the constant parameters for the calculation
    # surf_generator! is a function that generates the surface sp.ys for each iteration
    # N_ens is the number of ensemble averages to perform
    # returns the coherent and incoherent MDRC

    _, (pc_stats...) = @timed precompute!(data)
    @debug "Precompute stats: $pc_stats"

    thr_sz = Threads.nthreads()
    LinearAlgebra.BLAS.set_num_threads(thr_sz)

    sp = data.sp
    spa = data.spa
    iters = data.iters

    _, (mdrc_stats...) = @timed begin
    for n in ProgressBar(1:iters)
        generate_surface!(sp, spa)
        solve_single!(data)
        observe!(data, n)
    end
    end # end mdrc_stats

    @debug "MDRC stats: $mdrc_stats"

    return nothing
end

function get_coh_inc(out::SimOutput, qs, ks, above)
    coh = Matrix{Float64}(undef, (length(qs), length(ks)))
    inc = Matrix{Float64}(undef, (length(qs), length(ks)))
    for (j, k) in enumerate(ks)
        for (i, q) in enumerate(qs)
            C = abs2(alpha(q, above)) / abs(alpha(k, above))
            coh[i, j] = C * abs2(out.R[i, j])
            inc[i, j] = C * out.R²[i, j] - coh[i, j]
        end
    end
    return coh, inc
end

struct DataMDRC
    coh_p::Matrix{Float64}
    inc_p::Matrix{Float64}
    coh_s::Matrix{Float64}
    inc_s::Matrix{Float64}
end

"Returns coherent and incoherent MDRC calculated from ⟨R⟩ and ⟨R²⟩"
function get_qs_and_mdrc(data::SolverData)
    spa = data.spa
    ks = spa.ks
    qs = spa.qs

    reduced_qs = qs[qs .> -1.0 .&& qs .< 1.0]


    coh_p, inc_p = get_coh_inc(data.out_p, reduced_qs, ks, spa.above)
    coh_s, inc_s = get_coh_inc(data.out_s, reduced_qs, ks, spa.above)

    return reduced_qs, DataMDRC(coh_p, inc_p, coh_s, inc_s)
end
