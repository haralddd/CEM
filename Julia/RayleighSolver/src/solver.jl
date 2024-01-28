
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
    # Invariant under surface change and incident angle θ0

    ps = rp.ps
    qs = rp.qs
    κ = rp.ν == p ? rp.ε : rp.μ
    εμ = rp.ε * rp.μ

    @inbounds for n in axes(M, 3), j in axes(M, 2), i in axes(M, 1)
        p = ps[i]
        q = qs[j]
        M[i, j, n] = M_ker(p, q, κ, α(p, εμ), α0(q), n - 1)
    end
    return nothing
end

function N_invariant!(N::Array{ComplexF64,3}, rp::RayleighParams)::Nothing
    # Calculate the surface invariant part of the Npk vector
    # Invariant under surface change, but NOT incident angle θ0
    # Input ks must be the nearest value in qs for correct results
    ps = rp.ps
    qs = rp.qs
    κ = rp.ν == p ? rp.ε : rp.μ
    εμ = rp.ε * rp.μ

    @inbounds for n in axes(N, 3), (j, kj) in enumerate(rp.kis), i in axes(N, 1)
        p = ps[i]
        k = qs[kj]
        N[i, j, n] = N_ker(p, k, κ, α(p, εμ), α0(k), n - 1)
    end
    return nothing
end

function solve!(sp::SurfPreAlloc, rp::RayleighParams,
    M_pre::Array{ComplexF64,3}, N_pre::Array{ComplexF64,3})::Nothing

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

        for j in axes(sp.Mpq, 2), i in axes(sp.Mpq, 1)
            sp.Mpq[i, j] += M_pre[i, j, n] * sp.sFys[i+j-1]
        end

        for (j, kj) in enumerate(rp.kis), i in axes(sp.Npk, 1)
            sp.Npk[i, j] -= N_pre[i, j, n] * sp.sFys[i+kj-1]
        end
    end

    for j in eachindex(rp.kis)
        sp.Npk[:, j] .= sp.Mpq \ sp.Npk[:, j]
    end

    return nothing
end

function MDRC_prefactor(k, q, L)::Float64
    return 1 / (2π * L) * α0(q) / α0(k)
end

function solve_ensemble(rp::RayleighParams, sp::SurfPreAlloc, surf_generator!::Function, N_ens::Int)

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
    Npk_pre = Array{ComplexF64}(undef, length(rp.ps), length(rp.ks), rp.Ni + 1)
    @time N_invariant!(Npk_pre, rp)

    nan_idxs = findall(isnan.(Npk_pre))
    inf_idxs = findall(isinf.(Npk_pre))
    @assert length(nan_idxs) == 0 "Npk_pre has NaNs at indices $nan_idxs"
    @assert length(inf_idxs) == 0 "Npk_pre has Infs at indices $inf_idxs"

    # Reduced q-vector indices
    qis = findall(q -> q > -1 && q < 1, rp.qs)
    display("Solving for θ0, N: $N_ens")

    # Choose the surface type function
    coh = zeros(ComplexF64, length(qis))
    incoh = zeros(Float64, length(qis))

    @time for _ in 1:N_ens
        surf_generator!(sp.ys)
        solve!(sp, rp, Mpk_pre, Npk_pre)

        # Add to local variables
        coh += sp.Npk[qis]
        incoh += abs2.(sp.Npk[qis])
    end

    coh .= abs2.(coh ./ N_ens)
    incoh .= MDRC_prefactor.(k, rp.qs[qis], rp.Lx) .* (incoh ./ N_ens .- coh)
    coh .*= MDRC_prefactor.(k, k, rp.Lx)
    return coh, incoh
end