
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
    sFys_pkidxs = spa.sFys_pkidxs

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

        # display(findall(!iszero, sFys))
        # error("Debugging")

        for j in axes(pd.Mpq, 2)
            for i in axes(pd.Mpq, 1)
                idx = sFys_pqidxs[i, j]
                pd.Mpq[i, j] += pd.Mpqn[i, j, n] * sFys[idx]
                sd.Mpq[i, j] += sd.Mpqn[i, j, n] * sFys[idx]
            end
        end
        for j in axes(pd.Npk, 2)
            for i in axes(pd.Npk, 1)
                idx = sFys_pkidxs[i, j]
                pd.Npk[i, j] -= pd.Npkn[i, j, n] * sFys[idx]
                sd.Npk[i, j] -= sd.Npkn[i, j, n] * sFys[idx]
            end
        end
    end

    @debug findall(x->x != 0.0, pd.Mpq)
    # @debug pd.Mpq[findall(x -> x != 0.0, pd.Mpq)]

    A = lu!(pd.Mpq)
    for i in axes(pd.Npk, 2)
        b = @view pd.Npk[:, i]
        ldiv!(A, b)
    end

    A = lu!(sd.Mpq)
    for i in axes(sd.Npk, 2)
        b = @view sd.Npk[:, i]
        ldiv!(A, b)
    end

    # @info "Solve 1"
    # @time for i in axes(pd.Npk, 2)
    #     b = @view pd.Npk[:, i]
    #     b = pd.Mpq \ b
    # end

    # # @info "Solve 2"
    # @time for i in axes(sd.Npk, 2)
    #     b = @view sd.Npk[:, i]
    #     b = sd.Mpq \ b
    # end


    # @info "Factorize 1"
    # @time A = lu!(pd.Mpq)
    # @info typeof(A)

    # @info "Solve 1"
    # @time for i in axes(pd.Npk, 2)
    #     b = @view pd.Npk[:, i]
    #     ldiv!(A, b)
    # end

    # @info "Factorize 2"
    # @time A = lu!(sd.Mpq)
    # @info typeof(A)

    # @info "Solve 2"
    # @time for i in axes(sd.Npk, 2)
    #     b = @view sd.Npk[:, i]
    #     ldiv!(A, b)
    # end


    # @info "gesv!"
    # @time LinearAlgebra.LAPACK.gesv!(pd.Mpq, pd.Npk)
    # @time LinearAlgebra.LAPACK.gesv!(sd.Mpq, sd.Npk)

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

    out_p = data.out_p
    out_s = data.out_s

    out_p.R .= 0.0
    out_p.R² .= 0.0
    out_p.σ² .= 0.0
    out_p.κ .= 0.0
    out_s.R .= 0.0
    out_s.R² .= 0.0
    out_s.σ² .= 0.0
    out_s.κ .= 0.0

    @info "mainloop:"
    # for n in ProgressBar(1:iters)
    @time for n in 1:iters
        generate_surface!(sp, spa)
        solve_single!(data)
        # observe!(data, n)
        out_p.R .+= sp.p_data.Npk
        out_p.R² .+= abs2.(sp.p_data.Npk)
        out_s.R .+= sp.s_data.Npk
        out_s.R² .+= abs2.(sp.s_data.Npk)
    end

    out_p.R ./= iters
    out_p.R² ./= iters
    out_s.R ./= iters
    out_s.R² ./= iters

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

    return reduced_qs, DataMDRC(coh_p, inc_p, coh_s, inc_s)z
end
