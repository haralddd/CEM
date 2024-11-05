
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

struct SolverData
    spa::SimParams
    sp::SimPrealloc
    out::SimOutput
    iters::Int64

    function SolverData(spa::SimParams, iters::Int64)
        sp, (sp_stats...) = @timed SimPrealloc(spa, iters)
        out, (out_stats...) = @timed SimOutput(spa)

        @debug "SimOutput stats: $out_stats"
        @debug "SimPrealloc stats: $sp_stats"
        @debug "SimPreCompute stats: $pc_stats"

        return new(spa, pc, sp, out, iters)
    end
end

function precompute!(data::SolverData)::Nothing
    spa = data.spa
    pc = data.pc
    precompute!(pc, spa)
    validate(pc)
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
    pc = data.pc
    sp = data.sp

    Mpqn = pc.Mpqn
    Npkn = pc.Npkn

    FFT = spa.FFT

    rev_kis = spa.rev_kis
    rev_qis = spa.rev_qis

    Mpq = sp.Mpq
    Npk = sp.Npk
    ys = sp.ys
    Fys = sp.Fys
    sFys = sp.sFys


    Mpq .= 0.0
    Npk .= 0.0

    @inbounds for n in reverse(axes(Mpqn, 3)) # Reverse because prefactors vanish at higher powers of ´n´
        @inbounds for i in eachindex(Fys)
            Fys[i] = ys[i]^(n - 1)
        end

        FFT * Fys # In place FFT

        fftshift!(sFys, Fys)

        @inbounds for (j, qj) in enumerate(rev_qis)
            @inbounds for i in axes(Mpq, 1)
                Mpq[i, j] += Mpqn[i, j, n] * sFys[i+qj-1]
            end
        end

        @inbounds for (j, kj) in enumerate(rev_kis)
            @inbounds for i in axes(Npk, 1)
                Npk[i, j] -= Npkn[i, j, n] * sFys[i+kj-1]
            end
        end
    end

    LinearAlgebra.LAPACK.gesv!(Mpq, Npk)

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
    A = sp.Npk

    _, (mdrc_stats...) = @timed begin
    for n in ProgressBar(1:iters)
        generate_surface!(sp, spa)
        solve_single!(data)
        observe!(data.out, A, n)
    end
    end # end mdrc_stats

    @debug "MDRC stats: $mdrc_stats"

    return nothing
end
"Returns coherent and incoherent MDRC calculated from ⟨R⟩ and ⟨R²⟩"
function get_mdrc_qs_coh_inc(data::SolverData)::Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Float64}}
    out = data.out
    spa = data.spa
    ks = spa.ks
    qs = spa.qs

    R = out.R
    R² = out.R²

    reduced_qs = qs[qs .> -1.0 .&& qs .< 1.0]

    coh = Matrix{Float64}(undef, (length(reduced_qs), length(ks)))
    inc = Matrix{Float64}(undef, (length(reduced_qs), length(ks)))
    for (j, k) in enumerate(ks)
        for (i, q) in enumerate(reduced_qs)
            C = abs2(alpha(q, spa.above)) / abs(alpha(k, spa.above))
            coh[i, j] = C * abs2(R[i, j])
            inc[i, j] = C * R²[i, j] - coh[i, j]
        end
    end

    return reduced_qs, coh, inc
end
