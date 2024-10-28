
struct SimOutput
    qs::Vector{Float64}
    coherent::Matrix{ComplexF64}
    incoherent::Matrix{ComplexF64}

    function SimOutput(spa::SimParams)
        # Reduced q-vector indices
        qis = findall(q -> q > -1 && q < 1, spa.qs)
        new_qs = spa.qs[qis]

        coh = zeros(Float64, (length(new_qs), length(spa.kis)))
        incoh = zeros(Float64, (length(new_qs), length(spa.kis)))
        new(coh, incoh, new_qs)
    end
end

struct SolverData
    spa::SimParams
    pc::SimPreCompute
    sp::SimPrealloc
    out::SimOutput
    iters::Int64

    function SolverData(spa::SimParams, iters::Int64)
        sp = SimPrealloc(spa)
        out = SimOutput(spa)

        pc = SimPreCompute(spa)
        validate(pc)

        return new(spa, pc, sp, out, iters)
    end
end

"""
    function solve!(data::SolverData)::Nothing

Calculates the preallocated surface integral.
Matrix Mpqn is preallocated and contains the invariant parts of the Mpq matrix.
Matrix Npkn is preallocated and contains the invariant parts of the Npk vector.
Stores the resulting reflection coefficients in sp.Npk
Overwrites sp.Mpq with the LU factorization in the linear solution process.

# Arguments:
- `data`: [`SolverData`](@ref) - Contains the parameters, preallocated steps, output
"""
function solve!(data::SolverData)::Nothing

    spa = data.spa
    pc = data.pc
    sp = data.sp

    Mpqn = pc.Mpqn
    Npkn = pc.Npkn

    FFT = spa.FFT
    kis = spa.kis

    Mpq = sp.Mpq
    Npk = sp.Npk
    ys = sp.ys
    Fys = sp.Fys
    sFys = sp.sFys


    Mpq .= 0.0
    Npk .= 0.0
    # To access the the Fourier transform of the surface integral
    # we must access the pattern ζ(p-q), so make a reverse index range
    qidxs = reverse(axes(Mpq, 2))
    kidxs = reverse(kis)

    @inbounds for n in reverse(axes(Mpqn, 3)) # Reverse because prefactors vanish at higher powers of ´n´
        @inbounds for i in eachindex(Fys)
            Fys[i] = ys[i]^(n - 1)
        end

        FFT * Fys # In place FFT

        fftshift!(sFys, Fys)

        @inbounds for (j, qj) in enumerate(qidxs)
            @inbounds for i in axes(Mpq, 1)
                Mpq[i, j] += Mpqn[i, j, n] * sFys[i+qj-1]
            end
        end

        @inbounds for (j, kj) in enumerate(kidxs)
            @inbounds for i in axes(Npk, 1)
                Npk[i, j] -= Npkn[i, j, n] * sFys[i+kj-1]
            end
        end
    end

    Npk, Mpq, _ = LinearAlgebra.LAPACK.gesv!(Mpq, Npk)

    return nothing
end

function MDRC_prefactor(k::Float64, q::Float64, L::Float64)::Float64
    return 1 / (2π * L) * alpha0(q) / alpha0(k)
end

function solve_MDRC!(data::SolverData)
    # Solve the Rayleigh problem for a given SimParams struct
    # sp is the preallocated surface struct, containing the preallocated output, this is mutated in place
    # spa is the SimParams struct, containing the constant parameters for the calculation
    # surf_generator! is a function that generates the surface sp.ys for each iteration
    # N_ens is the number of ensemble averages to perform
    # returns the coherent and incoherent MDRC

    thr_sz = Threads.nthreads()
    LinearAlgebra.BLAS.set_num_threads(thr_sz)

    sp = data.sp
    spa = data.spa
    iters = data.iters
    qs_reduced = data.out.qs
    valid_qis = findall(q -> q > -1 && q < 1, spa.qs)
    coh = data.out.coherent
    incoh = data.out.incoherent

    @time begin
        for n in 1:iters
            generate_surface!(sp, spa)
            solve!(data)

            # Accumulate results
            for j in eachindex(spa.kis)
                for (i, qi) in enumerate(valid_qis)
                    coh[i, j] += sp.Npk[qi, j]
                    incoh[i, j] += abs2(sp.Mpq[qi, j])
                end
            end
        end
    end # time

    for (j, k) in enumerate(spa.ks)
        for (i, q) in enumerate(qs_reduced)
            prefactor = abs2(alpha0(q)) / abs(alpha0(k))
            incoh[i, j] = prefactor * mean(abs2.(res[i, j, :])) - coh[i, j]
            coh[i, j] = prefactor * abs2.(mean(res[i, j, :]))
        end
    end

    return qs_reduced, coh, incoh
end