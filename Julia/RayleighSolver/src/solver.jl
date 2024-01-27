
function α(q::Float64, εμ::ComplexF64)::ComplexF64
    # Calculate the α parameter for a given p
    return sqrt(εμ - q^2)
end

function α0(q::Float64)::ComplexF64
    # Calculate the α0 parameter for a given q
    return sqrt(complex(1.0 - q^2))
end

function M_ker(p::Float64, q::Float64, κ::ComplexF64, α::ComplexF64, α0::ComplexF64, n::Int)::ComplexF64
    # Calculate the kernel of the Mpq matrix
    Δα = α - α0
    return (-1.0im)^n / factorial(n) * (
        (p + κ * q) * (p - q) * Δα^(n - 1) +
        (α + κ * α0) * Δα^n
    )
end

function N_ker(p::Float64, k::Float64, κ::ComplexF64, α::ComplexF64, α0::ComplexF64, n::Int)::ComplexF64
    Δα = α + α0
    return (-1.0im)^n / factorial(n) * (
        (p + κ * k) * (p - k) * Δα^(n - 1) +
        (α - κ * α0) * Δα^n
    )
end

function M_invariant!(M::Array{ComplexF64,3}, rp::RayleighParams)::Nothing
    # Calculate the surface invariant part of the Mpq matrix
    # Invariant under surface change), but NOT incident angle θ0

    ps = rp.ps
    qs = rp.qs
    κ = rp.ν == p ? rp.ε : rp.μ
    εμ = rp.ε * rp.μ

    for n in axes(M, 3)
        for i in eachindex(ps)
            p = ps[i]
            for j in eachindex(qs)
                q = qs[j]
                M[i, j, n] = M_ker(p, q, κ, α(p, εμ), α0(q), n - 1)
            end
        end
    end
    return nothing
end

function N_invariant!(N::Matrix{ComplexF64}, rp::RayleighParams, k::Float64)::Nothing
    # Calculate the surface invariant part of the Npk vector
    # Invariant under surface change, but NOT incident angle θ0
    # Input k must be the nearest value in qs for correct results
    ps = rp.ps
    κ = rp.ν == p ? rp.ε : rp.μ
    εμ = rp.ε * rp.μ

    for n in axes(N, 2)
        for i in eachindex(ps)
            p = ps[i]
            N[i, n] = N_ker(p, k, κ, α(p, εμ), α0(k), n - 1)
        end
    end
    return nothing
end

function solve!(sp::SurfPreAlloc, rp::RayleighParams,
    M_pre::Array{ComplexF64,3}, N_pre::Matrix{ComplexF64},
    ki::Int)::Nothing

    # Calculates the preallocated surface integral
    # sp is the preallocated surface struct, containing the preallocated output
    # rp is the RayleighParams struct, containing the constant parameters for the calculation
    # Matrix M_pre is preallocated and contains the invariant part of the Mpq matrix
    # Matrix N_pre is preallocated and contains the invariant part of the Npk vector
    # ki is the index of the nearest value in qs to the incident angle θ0

    # Stores the result in sp.Npk

    sp.Mpq .= 0.0
    sp.Npk .= 0.0

    for n in axes(M_pre, 3)
        sp.Fys .= sp.ys .^ (n - 1)
        rp.FT * sp.Fys # In place FFT
        fftshift!(sp.sFys, sp.Fys)

        for i in axes(sp.Mpq, 1), j in axes(sp.Mpq, 2)
            sp.Mpq[i, j] += M_pre[i, j, n] * sp.sFys[i+j-1]
        end

        for i in eachindex(sp.Npk)
            sp.Npk[i] -= N_pre[i, n] * sp.sFys[i+ki-1]
        end
    end

    sp.Npk .= sp.Mpq \ sp.Npk
    return nothing
end

function MDRC_prefactor(k, q, L)::Float64
    return 1 / (2π * L) * α0(q) / α0(k)
end

function run_threaded(rp; θ=10.0, N_ens::Int=10, surf_t::SurfType=gaussian)

    # Calc the invariant part of Mpk
    display("Calculating invariant parts of Mpk")
    Mpk_pre = Array{ComplexF64,3}(undef, length(rp.ps), length(rp.qs), rp.Ni + 1)
    @time M_invariant!(Mpk_pre, rp)

    # Check for NaNs and Infs
    nan_idxs = findall(isnan.(Mpk_pre))
    inf_idxs = findall(isinf.(Mpk_pre))
    @assert length(nan_idxs) == 0 "Mpk_pre has NaNs at indices $nan_idxs"
    @assert length(inf_idxs) == 0 "Mpk_pre has Infs at indices $inf_idxs"

    #### First angle
    ki = searchsortedfirst(rp.qs, sind(θ), rev=true)
    k = rp.qs[ki]

    display("Calculating invariant parts of Npk")
    Npk_pre = Matrix{ComplexF64}(undef, length(rp.ps), rp.Ni + 1)
    @time N_invariant!(Npk_pre, rp, k)

    nan_idxs = findall(isnan.(Npk_pre))
    inf_idxs = findall(isinf.(Npk_pre))
    @assert length(nan_idxs) == 0 "Npk_pre has NaNs at indices $nan_idxs"
    @assert length(inf_idxs) == 0 "Npk_pre has Infs at indices $inf_idxs"

    # Reduced q-vector indices
    qis = findall(q -> q > -1 && q < 1, rp.qs)
    display("Solving for θ0, N: $N_ens")

    # Vector of atomic variables to accumulate results
    coh_re = [Atomic{Float64}(0.0) for _ in qis]
    coh_im = [Atomic{Float64}(0.0) for _ in qis]
    inc = [Atomic{Float64}(0.0) for _ in qis]

    # Write safe preallocation of surface structs
    T = nthreads()
    sps = [SurfPreAlloc(rp, surf_t) for _ in 1:T]

    @time @threads for t in 1:T
        # Get local variables
        sp = sps[t] # Surface prealloc struct is mutated in place

        coh_local = Vector{ComplexF64}(undef, length(qis)) # Accumulate coherent part
        inc_local = Vector{Float64}(undef, length(qis)) # Accumulate term to subtract from coherent part

        blk = N_ens ÷ T + (t <= N_ens % T ? 1 : 0) # Number of realizations to solve for in this thread

        for _ in 1:blk # Distribute workload of N_ens over threads
            generate!(sp.ys, rp, surf_t)
            solve!(sp, rp, Mpk_pre, Npk_pre, ki)

            # Add to local variables
            coh_local .+= sp.Npk[qis]
            inc_local .+= abs2.(sp.Npk[qis]) # First accumulate the term to subtract from the coherent part
        end

        # Critical section could be better here
        for i in eachindex(coh_re)
            atomic_add!(coh_re[i], real(coh_local[i]))
        end
        for i in eachindex(coh_im)
            atomic_add!(coh_im[i], imag(coh_local[i]))
        end
        for i in eachindex(inc)
            atomic_add!(inc[i], inc_local[i])
        end
    end

    # Take mean and find incoherent part
    coh = Vector{Float64}(undef, length(qis))
    incoh = Vector{Float64}(undef, length(qis))

    for (i, qi) in enumerate(qis)

        coh[i] = ((coh_re[i][])^2 + (coh_im[i][])^2) / N_ens
        incoh[i] = coh[i] - (inc[i][] / N_ens)

        pre = MDRC_prefactor(k, rp.qs[qi], L)

        coh[i] *= pre
        incoh[i] *= pre

    end
    return coh, incoh
end