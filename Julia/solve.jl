# using MKL
push!(LOAD_PATH, "Julia/RayleighSolver/")
using RayleighSolver
using DelimitedFiles
using Dates
using Statistics

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

function MDRC_prefactor(k, q, L)
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

function solve_gaussian_and_save(; ν::Polarization=s, ε=2.25, μ=1.0, N_ens::Int=10)

    λ = 632.8e-9 # He-Ne laser wavelength
    L = 50 * λ # Length of surface
    δ = λ / 20
    a = λ / 4

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
    end

    # Ensemble params
    θ0 = 0.0
    θ1 = 34.05

    @assert N_ens != Nq "N_ens must be different from Nq due to binary file read limitations"

    # Calc the invariant part of Mpk
    display("Calculating invariant parts of Mpk")
    Mpk_pre = Array{ComplexF64,3}(undef, length(rp.ps), length(rp.qs), Ni + 1)
    @time pre_M_invariant!(Mpk_pre, rp)
    nan_idxs = findall(isnan.(Mpk_pre))
    inf_idxs = findall(isinf.(Mpk_pre))
    @assert length(nan_idxs) == 0 "Mpk_pre has NaNs at indices $nan_idxs"
    @assert length(inf_idxs) == 0 "Mpk_pre has Infs at indices $inf_idxs"


    # display(Mpk_pre[isnan.(Mpk_pre)])
    # display(Mpk_pre[isinf.(Mpk_pre)])

    res = Matrix{ComplexF64}(undef, length(rp.qs), N_ens)

    #### First angle
    ki = searchsortedfirst(rp.qs, sind(θ0), rev=true)
    k = rp.qs[ki]

    display("Calculating invariant parts of Npk")
    Npk_pre = Matrix{ComplexF64}(undef, length(rp.ps), Ni + 1)
    @time pre_N_invariant!(Npk_pre, rp, k)

    nan_idxs = findall(isnan.(Npk_pre))
    inf_idxs = findall(isinf.(Npk_pre))
    @assert length(nan_idxs) == 0 "Npk_pre has NaNs at indices $nan_idxs"
    @assert length(inf_idxs) == 0 "Npk_pre has Infs at indices $inf_idxs"

    display("Solving for θ0, N: $N_ens")
    sp = SurfPreAlloc(rp, gaussian)
    @time for i in axes(res, 2)
        sp = SurfPreAlloc(rp, gaussian)
        solve_pre!(sp, rp, Mpk_pre, Npk_pre, ki)
        res[:, i] .= sp.R
    end

    filestr = run_dir * "/θ$(θ0)_R_ComplexF64.bin"
    open(filestr, "w") do io
        write(io, res)
        display("Wrote to $filestr")
    end # io θ0

    # Save MDRC_incoh and MDRC_coh
    open(run_dir * "/θ$(θ0)_mdrc_incoh.bin", "w") do io
        mdrc_incoh = MDRC_incoh(res, k, rp.qs, L)
        write(io, mdrc_incoh)
    end
    open(run_dir * "/θ$(θ0)_mdrc_coh.bin", "w") do io
        mdrc_coh = MDRC_coh(res, k, rp.qs, L)
        write(io, mdrc_coh)
    end

    #### Second angle

    ki = searchsortedfirst(rp.qs, sind(θ1), rev=true)
    k = rp.qs[ki]

    display("Calculating invariant parts of Npk")
    pre_N_invariant!(Npk_pre, rp, k)
    display("Solving for θ1, N: $N_ens")
    @time for i in axes(res, 2)
        sp = SurfPreAlloc(rp, gaussian)
        solve_pre!(sp, rp, Mpk_pre, Npk_pre, ki)
        res[:, i] .= sp.R
    end

    filestr = run_dir * "/θ$(θ1)_R_ComplexF64.bin"
    open(filestr, "w") do io
        write(io, res)
        display("Wrote to $filestr")
    end # io θ1

    # Save MDRC_incoh and MDRC_coh
    open(run_dir * "/θ$(θ1)_mdrc_incoh.bin", "w") do io
        mdrc_incoh = MDRC_incoh(res, k, rp.qs, L)
        write(io, mdrc_incoh)
    end
    open(run_dir * "/θ$(θ1)_mdrc_coh.bin", "w") do io
        mdrc_coh = MDRC_coh(res, k, rp.qs, L)
        write(io, mdrc_coh)
    end
    nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 4
        println("Usage: julia solve.jl [p|s] [ε] [μ] [N].
        p|s:\tPolarization
        ε:\tRelative permittivity of the material
        μ:\tRelative permeability of the material
        N:\tNumber of surface realizations to solve for")
        exit(1)
    end
    pol = polarization_from_string(ARGS[1])
    ε = parse(ComplexF64, ARGS[2])
    μ = parse(ComplexF64, ARGS[3])
    N = parse(Int, ARGS[4])

    solve_gaussian_and_save(; ν=pol, ε=ε, μ=μ, N_ens=N)
else

end