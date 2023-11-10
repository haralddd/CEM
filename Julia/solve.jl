push!(LOAD_PATH, "Julia/RayleighSetup/")

display(Base.load_path())
using RayleighSetup


function solve_fresnel()
    surf_t::SurfType = flat
    RayleighSetup.rp(
        ν=p,
        Nq=2^11,
        ε=2.25,
        L=10.0e-6,
        Q_mult=4,
        Ni=10
    )


end


function solve_surf(rp::RayleighParams, sp::SurfPreAlloc)

    for n in 0:rp.Ni
        sp.Fys .= ys .^ n
        rp.FT, rp.Fys # In place FFT
        fftshift!(rp.sFys, rp.Fys) # No way to do shift in place apparently


        for I in CartesianIndices(rp.Mpq)
            i, j = Tuple(I)
            rp.Mpq[I] += M_pre[i, j, n] * rp.sFys[i+j-1]
        end

        for i in eachindex(rp.Npk)
            rp.Npk[i] -= N_pre[i, n] * rp.sFys[i+rp.ki-1]
        end
    end

    return rp.Mpq \ rp.Npk
end



function solve(rp::RayleighParams, reduction::Function=identity)
    # Solve the reduced Rayleigh equations for p-polarized light

    # Shorthand parameters
    κ = rp.ν == p ? rp.ε : rp.μ
    k = rp.k
    ki = rp.ki
    Ni = rp.Ni

    α(q::Float64)::ComplexF64 = √complex(rp.ε * rp.μ - q^2)

    # Assumed μ0 = ε0 = 1
    α0(q::Float64)::ComplexF64 = √complex(1.0 - q^2)

    Mpq_eval(p::Float64, q::Float64, n::Int)::ComplexF64 = (
        (-1.0im)^n / factorial(n) * (
            (p + κ * q) * (p - q) * (α(p) - α0(q))^(n - 1) +
            (α(p) + κ * α0(q)) * (α(p) - α0(q))^n)
    )
    Npk_eval(p::Float64, n::Int)::ComplexF64 = (
        (-1.0im)^n / factorial(n) * (
            (p + κ * k) * (p - k) * (α(p) + α0(k))^(n - 1) +
            (α(p) - κ * α0(k)) * (α(p) + α0(k))^n)
    )

    for n in 0:Ni
        display("n = $n")

        rp.Fys .= rp.ys .^ n
        rp.FT * rp.Fys # In place FFT
        fftshift!(rp.sFys, rp.Fys) # No way to do shift in place apparently

        for I in CartesianIndices(rp.Mpq)
            i, j = Tuple(I)
            rp.Mpq[I] += Mpq_eval(rp.ps[i], rp.qs[j], n) * rp.sFys[i+j-1]
        end

        for i in eachindex(rp.Npk)
            rp.Npk[i] -= Npk_eval(rp.ps[i], n) * rp.sFys[i+ki-1]
        end
    end

    nans1 = findall(!isfinite, rp.Mpq)
    nans2 = findall(!isfinite, rp.Npk)

    if length(nans1) > 0 || length(nans2) > 0
        println("Mpq:")
        display(rp.Mpq[nans1])
        println("Mpq idxs:")
        display(nans1)

        println("Npk:")
        display(rp.Npk[nans2])
        println("Npk idxs:")
        display(nans2)

        display("Number of NaNs in Mpq: $(length(nans1))")
        display("Number of NaNs in Npk: $(length(nans2))")

        throw("Non-finite values in Mpq or Npk")
    end

    # Solve the system Δq / 2π ∑_q N⁺ₚ(p|q) * Rₚ(q|k) = N⁻ₚ(p|k) for R
    # rp.Rqk .= rp.Mpq \ rp.Npk

    return reduction(rp.Mpq \ rp.Npk)
end

mean_DRC(rp::RayleighParams, sol::Matrix{ComplexF64}) = begin
    # Calculate the mean differential reflection coefficient
    # Factor of ω / c0 removed due to scaling

    ω = rp.ω
    c0 = rp.c0
    qs = rp.qs[rp.qs.<ω/c0.&&rp.qs.>-ω/c0]
    ks = rp.ks

    retval = Matrix{Float64}(undef, length(qs), length(ks))

    for (i, q) in enumerate(qs)
        for (j, k) in enumerate(ks)
            retval[i, j] = 1.0 / rp.Lx / 2π * cos(θ(q, ω, c0))^2 / cos(θ(k, ω, c0)) * abs2(sol[i, j])
        end
    end
    return retval
end