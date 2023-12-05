using MKL
push!(LOAD_PATH, "Julia/RayleighSetup/")
push!(LOAD_PATH, "Julia/RayleighSolver/")
using RayleighSetup
using RayleighSolver
using DelimitedFiles
using Dates

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

ε_list = Dict("glass" => 2.25 + 1.0e-4im, "silver" => -17.5 + 0.48im)
function solve_gaussian(; type="glass", ν::Polarization=s, N_ens::Int=10)

    λ = 632.8e-9 # He-Ne laser wavelength
    L = 50 * λ # Length of surface
    δ = λ / 20
    a = λ / 4

    Nq = 2^10

    ε = ε_list[type]
    Ni = 10


    ### Solve for p-polarization
    # General parameters
    rp = RayleighParams(
        ν=ν,
        λ=λ,
        Nq=Nq,
        ε=ε,
        Ni=Ni,
        L=L,
        δ=δ,
        a=a,
    )
    display(rp.ν)


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
    @time for i in axes(res, 2)
        sp = SurfPreAlloc(rp, gaussian)
        solve_pre!(sp, rp, Mpk_pre, Npk_pre, ki)
        res[:, i] .= sp.R
    end

    timestamp = now() |> string
    run_dir = "data/" * timestamp
    setup_dir(run_dir)
    filestr = run_dir * "/gsilver_θ$(θ0)_ν$(ν|>string)_Nq$(Nq)_Nens$(N_ens).bin"
    open(filestr, "w") do io
        write(io, res)
        display("Wrote to $filestr")
    end # io θ0

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

    filestr = run_dir * "/gsilver_θ$(θ1)_ν$(ν|>string)_Nq$(Nq)_Nens$(N_ens).bin"
    open(filestr, "w") do io
        write(io, res)
        display("Wrote to $filestr")
    end # io θ1
    nothing
end

if length(ARGS) != 3
    println("Usage: julia solve.jl [glass|silver] [p|s] [N] . Where p|s is the polarization and N is the number of surfaces to solve for.")
    exit(1)
end

type = ARGS[1]
pol = polarization_from_string(ARGS[2])
N = parse(Int, ARGS[3])

solve_gaussian(; type=type, ν=pol, N_ens=N)
exit(0)