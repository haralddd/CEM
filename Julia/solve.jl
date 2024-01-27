# using MKL
push!(LOAD_PATH, "Julia/RayleighSolver/")
using RayleighSolver
using Dates
using Statistics
using LinearAlgebra
using Base.Threads
LinearAlgebra.BLAS.set_num_threads(1)

function setup_dir(str="data")
    if !isdir(str)
        mkdir(str)
    end
end

function test_write()
    setup_dir()
    timestamp = now() |> string
    @show A1 = [1.0 + 1.0im, 2.0 + 2.0im, 0] .|> complex
    @show A2 = [3.0 + 3.0im, 3, 4.0im,] .|> complex
    N = 3
    pol = s
    open("data/" * timestamp * "test_$(pol|>string)_N$(N).bin", "w") do io
        write(io, A1)
        write(io, A2)
    end
    open("data/" * timestamp * "test_$(pol|>string)_N$(N).bin", "r") do io
        data = reinterpret(ComplexF64, read(io))
        data = reshape(data, 3, :)
        display(data)
        display(data[:, 1])
    end
end

function MDRC_prefactor(k, q, L)::Float64
    return 1 / (2π * L) * α0(q) / α0(k)
end

function MDRC_coh(R, k, qs, L)
    # Find indices of qs where qs is between -1 and 1
    qis = findall(q -> q > -1 && q < 1, qs)
    mdrc_coh = Vector{Float64}(undef, length(qis))

    for (i, qi) in enumerate(qis)
        mdrc_coh[i] = MDRC_prefactor(k, qs[qi], L) * abs2(mean(R[qi, :]))
    end
    return mdrc_coh
end

function MDRC_incoh(R, k, qs, L)
    # Find indices of qs where qs is between -1 and 1
    qis = findall(q -> q > -1 && q < 1, qs)
    mdrc_incoh = Vector{Float64}(undef, length(qis))
    for (i, qi) in enumerate(qis)
        mdrc_incoh[i] = MDRC_prefactor(k, qs[qi], L) * (mean(abs2.(R[qi, :])) - abs2(mean(R[qi, :])))
    end
    return mdrc_incoh
end

function solve_gaussian_and_save(; ν::Polarization=s, ε=2.25, μ=1.0, θ=10.0, N_ens::Int=10)

    λ = 457.9e-9 # Wavelength of incident light
    L = 100 * λ # Length of surface
    δ = 5.0e-9 # RMS height of surface
    a = 100e-9 # Correlation length of surface

    Nq = 2^10
    Ni = 10


    ### Solve for p-polarization
    # General parameters
    rp = RayleighParams(
        ν=ν,
        λ=λ,
        Nq=Nq,
        ε=ε,
        μ=μ,
        Ni=Ni,
        L=L,
        δ=δ,
        a=a,
    )

    # Directory setup and dump config to txt
    timestamp = now() |> string
    run_dir = "data/" * timestamp
    setup_dir(run_dir)
    open(run_dir * "/config.txt", "w") do io
        write(io, params_as_string(rp))
        write(io, "\nN_ens = $N_ens")
        write(io, "\nθ = $θ")
    end

    # Calc the invariant part of Mpk
    display("Calculating invariant parts of Mpk")
    Mpk_pre = Array{ComplexF64,3}(undef, length(rp.ps), length(rp.qs), Ni + 1)
    @time M_invariant!(Mpk_pre, rp)
    nan_idxs = findall(isnan.(Mpk_pre))
    inf_idxs = findall(isinf.(Mpk_pre))
    @assert length(nan_idxs) == 0 "Mpk_pre has NaNs at indices $nan_idxs"
    @assert length(inf_idxs) == 0 "Mpk_pre has Infs at indices $inf_idxs"

    #### First angle
    ki = searchsortedfirst(rp.qs, sind(θ), rev=true)
    k = rp.qs[ki]

    display("Calculating invariant parts of Npk")
    Npk_pre = Matrix{ComplexF64}(undef, length(rp.ps), Ni + 1)
    @time N_invariant!(Npk_pre, rp, k)
    sps = [SurfPreAlloc(rp, gaussian) for _ in 1:Threads.nthreads()]

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

    T = nthreads()

    @threads for t in 1:T
        # Get local variables
        sp = sps[t] # Surface prealloc struct is mutated in place

        coh_local = Vector{ComplexF64}(undef, length(qis)) # Accumulate coherent part
        inc_local = Vector{Float64}(undef, length(qis)) # Accumulate term to subtract from coherent part

        blk = N_ens ÷ T + (t <= N_ens % T ? 1 : 0) # Number of realizations to solve for in this thread

        for _ in 1:blk # Distribute workload of N_ens over threads
            generate!(sp.ys, rp, gaussian)
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

    # Save MDRC_incoh and MDRC_coh
    open(run_dir * "/incoh.bin", "w") do io
        write(io, incoh)
    end
    open(run_dir * "/coh.bin", "w") do io
        write(io, coh)
    end
    nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 5
        println("Usage: julia solve.jl [p|s] [ε] [μ] [N] [θ].
        p|s:\tPolarization
        ε:\tRelative permittivity of the material
        μ:\tRelative permeability of the material
        N:\tNumber of surface realizations to solve for
        θ:\tIncident angle in degrees")
        exit(1)
    end
    pol = polarization_from_string(ARGS[1])
    ε = parse(ComplexF64, ARGS[2])
    μ = parse(ComplexF64, ARGS[3])
    N = parse(Int, ARGS[4])
    θ = parse(Float64, ARGS[5])

    solve_gaussian_and_save(; ν=pol, ε=ε, μ=μ, N_ens=N, θ=θ)
else
    # Mostly for testing
    pol = s
    ε = -7.5 + 0.24im
    μ = 1.0
    N = 10000
    θ = 20.0
    solve_gaussian_and_save(; ν=pol, ε=ε, μ=μ, N_ens=N, θ=θ)
end