



"""
    function solve!(sp::SimPrealloc, spa::SimParams, pc::SimPreCompute)::Nothing

Calculates the preallocated surface integral.
Matrix M_pre is preallocated and contains the invariant part of the Mpq matrix.
Matrix N_pre is preallocated and contains the invariant part of the Npk vector.
ki is the index of the nearest value in qs to the incident angle θ0.
Stores the result in sp.Npk.

# Arguments:
- `sp`: [`SimPrealloc`](@ref) - Contains the preallocated steps and output, all fields are mutated
- `spa`: [`SimParams`](@ref) - Contains the constant parameters for the calculation
"""
function solve!(sp::SimPrealloc, spa::SimParams, pc::SimPreCompute)::Nothing


    Mpqn = pc.Mpqn
    Npkn = pc.Npkn

    Mpq = sp.Mpq
    Npk = sp.Npk

    FFT = spa.FFT
    ys = sp.ys
    Fys = sp.Fys
    sFys = sp.sFys

    kis = spa.kis

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

        fftshift!(sp.sFys, sp.Fys)

        @inbounds for (j, qj) in enumerate(qidxs)
            @inbounds for i in axes(sp.Mpq, 1)
                Mpq[i, j] += Mpqn[i, j, n] * sFys[i+qj-1]
            end
        end

        @inbounds for (j, kj) in enumerate(kidxs)
            @inbounds for i in axes(Npk, 1)
                Npk[i, j] -= Npkn[i, j, n] * sFys[i+kj-1]
            end
        end
    end

    LinearAlgebra.LAPACK.gesv!(Mpq, Npk)

    return nothing
end

function MDRC_prefactor(k::Float64, q::Float64, L::Float64)::Float64
    return 1 / (2π * L) * alpha0(q) / alpha0(k)
end



function solve_MDRC!(sp::SimPrealloc, spa::SimParams, N_ens::Int)
    # Solve the Rayleigh problem for a given SimParams struct
    # sp is the preallocated surface struct, containing the preallocated output, this is mutated in place
    # spa is the SimParams struct, containing the constant parameters for the calculation
    # surf_generator! is a function that generates the surface sp.ys for each iteration
    # N_ens is the number of ensemble averages to perform
    # returns the coherent and incoherent MDRC

    # Calc the invariant part of Mpk
    pc = SimPreCompute(spa)
    validate(pc)

    # Reduced q-vector indices
    qis = findall(q -> q > -1 && q < 1, spa.qs)
    qs_reduced = spa.qs[qis]
    thr_sz = Threads.nthreads()
    LinearAlgebra.BLAS.set_num_threads(thr_sz)

    # sp_vec = [deepcopy(sp) for _ in 1:thr_sz]
    res = zeros(ComplexF64, length(qis), length(spa.kis), N_ens)
    @time begin
        # Threads.@threads for n in 1:N_ens
        #     tid = Threads.threadid()
        #     my_sp = sp_vec[tid]
        #     generate!(my_sp, spa)
        #     solve!(my_sp, spa, Mpk_pre, Npk_pre)
        for n in 1:N_ens

            generate_surface!(sp, spa)
            solve!(sp, spa, pc)

            # Copy to local results
            for j in eachindex(spa.kis)
                for (i, qi) in enumerate(qis)
                    res[i, j, n] = sp.Npk[qi, j]
                end
            end
        end
    end # time

    coh = zeros(Float64, (length(qis), length(spa.kis)))
    incoh = zeros(Float64, (length(qis), length(spa.kis)))

    for (j, k) in enumerate(spa.ks)
        for (i, q) in enumerate(qs_reduced)
            prefactor = abs2(alpha0(q)) / abs(alpha0(k))
            coh[i, j] = prefactor * abs2.(mean(res[i, j, :]))
            incoh[i, j] = prefactor * mean(abs2.(res[i, j, :])) - coh[i, j]
        end
    end

    return qs_reduced, coh, incoh
end