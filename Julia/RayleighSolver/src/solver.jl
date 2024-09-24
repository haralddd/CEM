
"""
    function alpha(q::Float64, epsmu::ComplexF64)::ComplexF64

Calculate ``\\alpha = \\sqrt{\\epsilon \\mu - q^2}``

# Arguments:
- `q`: Wavenumber - Dimensionless, scaled by ``\\frac{\\omega}{c}``
- `epsmu`: ``\\epsilon\\cdot\\mu`` - pre-computed product of relative Permittivity and Permeability
"""
function alpha(q::Float64, epsmu::ComplexF64)::ComplexF64
    return sqrt(epsmu - q^2)
end

"""
    function alpha0(q::Float64)::ComplexF64

Calculate ``\\alpha_0 = \\sqrt{\\epsilon_0 \\mu_0 - q^2}``, where ``\\epsilon_0 = \\mu_0 = 1``

# Arguments:
- `q`: Wavenumber - Dimensionless, scaled by ``\\frac{\\omega}{c}``
"""
function alpha0(q::Float64)::ComplexF64
    return sqrt(complex(1.0 - q^2))
end

"""
    function M_ker(p::Float64, q::Float64, kappa::ComplexF64, alpha::ComplexF64, alpha0::ComplexF64, n::Int)::ComplexF64

Calculate the kernel of the Mpq matrix

# Arguments:
- p: Dimensionless, scaled by ``\\frac{\\omega}{c}``
- q: Dimensionless, scaled by ``\\frac{\\omega}{c}``
- kappa: Either ``\\varepsilon`` or ``\\mu``
- alpha: Wavenumber output from [`alpha`](@ref), scaled by ``\\frac{\\omega}{c}``
- alpha0: Wavenumber output from [`alpha0`](@ref), scaled by ``\\frac{\\omega}{c}``
- n: Current order in the Taylor expansion of the integral term ``I(\\gamma | q) = \\int_{-\\infty}^{\\infty} e^{-i\\gamma \\zeta x} e^{-iqx_1} dx``
"""
function M_ker(p::Float64, q::Float64, kappa::ComplexF64, alpha::ComplexF64, alpha0::ComplexF64, n::Int)::ComplexF64
    da = alpha - alpha0
    return (-1.0im)^n / factorial(n) * (
        (p + kappa * q) * (p - q) * da^(n - 1) +
        (alpha + kappa * alpha0) * da^n
    )
end


"""
    function N_ker(p::Float64, k::Float64, kappa::ComplexF64, alpha::ComplexF64, alpha0::ComplexF64, n::Int)::ComplexF64

Calculate the kernel of the Npk vector, see [`M_ker`](@ref), as the functionality is nearly identical.
"""
function N_ker(p::Float64, k::Float64, kappa::ComplexF64, alpha::ComplexF64, alpha0::ComplexF64, n::Int)::ComplexF64
    da = alpha + alpha0
    return (-1.0im)^n / factorial(n) * (
        (p + kappa * k) * (p - k) * da^(n - 1) +
        (alpha - kappa * alpha0) * da^n
    )
end

function M_invariant!(M::Array{ComplexF64,3}, rp::RayleighParams)::Nothing
    # Calculate the surface invariant part of the Mpq matrix
    # Invariant under surface change and incident angle θ0

    ps = rp.ps
    qs = rp.qs
    kappa = rp.nu == p ? rp.eps : rp.mu
    epsmu = rp.eps * rp.mu

    @inbounds for n in axes(M, 3), j in axes(M, 2), i in axes(M, 1)
        p = ps[i]
        q = qs[j]
        M[i, j, n] = M_ker(p, q, kappa, alpha(p, epsmu), alpha0(q), n - 1)
    end
    return nothing
end

function N_invariant!(N::Array{ComplexF64,3}, rp::RayleighParams)::Nothing
    # Calculate the surface invariant part of the Npk vector
    # Invariant under surface change, but NOT incident angle θ0
    # Input ks must be the nearest value in qs for correct results
    ps = rp.ps
    qs = rp.qs
    kappa = rp.nu == p ? rp.eps : rp.mu
    epsmu = rp.eps * rp.mu

    @inbounds for n in axes(N, 3), (j, kj) in enumerate(rp.kis), i in axes(N, 1)
        p = ps[i]
        k = qs[kj]
        N[i, j, n] = N_ker(p, k, kappa, alpha(p, epsmu), alpha0(k), n - 1)
    end
    return nothing
end

"""
    function solve!(sp::SimulationPreAlloc, rp::RayleighParams)

Calculates the preallocated surface integral.
Matrix M_pre is preallocated and contains the invariant part of the Mpq matrix.
Matrix N_pre is preallocated and contains the invariant part of the Npk vector.
ki is the index of the nearest value in qs to the incident angle θ0.
Stores the result in sp.Npk.

# Arguments:
- `sp`: [`SimulationPreAlloc`](@ref) - Contains the preallocated steps and output, all fields are mutated
- `rp`: [`RayleighParams`](@ref) - Contains the constant parameters for the calculation
"""
function solve!(sp::SimulationPreAlloc, rp::RayleighParams,
    M_pre::Array{ComplexF64,3}, N_pre::Array{ComplexF64,3})::Nothing


    sp.Mpq .= 0.0
    sp.Npk .= 0.0
    # To access the the Fourier transform of the surface integral
    # we must access the pattern ζ(p-q), so make a reverse index range
    qidxs = reverse(axes(sp.Mpq, 2))
    kidxs = reverse(rp.kis)

    for n in reverse(axes(M_pre, 3)) # Reverse because prefactors vanish at higher powers of ´n´
        for i in eachindex(sp.Fys)
            sp.Fys[i] = sp.ys[i] ^ (n - 1)
        end

        rp.FFT * sp.Fys # In place FFT

        fftshift!(sp.sFys, sp.Fys)

        for (j, qj) in enumerate(qidxs)
            for i in axes(sp.Mpq, 1)
                sp.Mpq[i, j] += M_pre[i, j, n] * sp.sFys[i+qj-1]
            end
        end

        for (j, kj) in enumerate(kidxs)
            for i in axes(sp.Npk, 1)
                sp.Npk[i, j] -= N_pre[i, j, n] * sp.sFys[i+kj-1]
            end
        end
    end

    LinearAlgebra.LAPACK.gesv!(sp.Mpq, sp.Npk)

    return nothing
end

function MDRC_prefactor(k::Float64, q::Float64, L::Float64)::Float64
    return 1 / (2π * L) * alpha0(q) / alpha0(k)
end

function precalc(rp)
    Mpk_pre = Array{ComplexF64,3}(undef, length(rp.ps), length(rp.qs), rp.Ni + 1)
    Npk_pre = Array{ComplexF64,3}(undef, length(rp.ps), length(rp.kis), rp.Ni + 1)

    M_invariant!(Mpk_pre, rp)
    N_invariant!(Npk_pre, rp)

    return Mpk_pre, Npk_pre
end

function solve_MDRC!(rp::RayleighParams, sp::SimulationPreAlloc, N_ens::Int)
    # Solve the Rayleigh problem for a given RayleighParams struct
    # sp is the preallocated surface struct, containing the preallocated output, this is mutated in place
    # rp is the RayleighParams struct, containing the constant parameters for the calculation
    # surf_generator! is a function that generates the surface sp.ys for each iteration
    # N_ens is the number of ensemble averages to perform
    # returns the coherent and incoherent MDRC

    # Calc the invariant part of Mpk
    @info "Calculating invariant parts of Mpk"
    @time Mpk_pre, Npk_pre = precalc(rp)

    @assert all(isfinite.(Mpk_pre))
    @assert all(isfinite.(Npk_pre))

    # Reduced q-vector indices
    qis = findall(q -> q > -1 && q < 1, rp.qs)
    qs_reduced = rp.qs[qis]
    @info "Solving for R, N: $N_ens"


    @show thr_sz = Threads.nthreads()
    LinearAlgebra.BLAS.set_num_threads(thr_sz)

    sp_vec = [deepcopy(sp) for _ in 1:thr_sz]
    res = zeros(ComplexF64, length(qis), length(rp.kis), N_ens)
    @time begin
    # Threads.@threads for n in 1:N_ens
    #     tid = Threads.threadid()
    #     my_sp = sp_vec[tid]
    #     generate!(my_sp, rp)
    #     solve!(my_sp, rp, Mpk_pre, Npk_pre)
    for n in 1:N_ens

        generate!(sp, rp)
        solve!(sp, rp, Mpk_pre, Npk_pre)

        # Add to local variables
        for j in eachindex(rp.kis)
            for (i, qi) in enumerate(qis)
                res[i, j, n] = sp.Npk[qi, j]
            end
        end
    end
    end # time

    coh = zeros(Float64, (length(qis), length(rp.kis)))
    incoh = zeros(Float64, (length(qis), length(rp.kis)))

    for (j, k) in enumerate(rp.ks)
        for (i, q) in enumerate(qs_reduced)
            prefactor = abs2(alpha0(q)) / abs(alpha0(k))
            coh[i, j] = prefactor * abs2.(mean(res[i, j, :]))
            incoh[i, j] = prefactor * mean(abs2.(res[i, j, :])) - coh[i, j]
        end
    end

    return qs_reduced, coh, incoh
end