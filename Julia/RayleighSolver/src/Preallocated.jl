# Preallocated structs which are mutated in place. These should be thread local.

struct Results
    R::Matrix{ComplexF64}
    R²::Matrix{Float64}
    σ²::Matrix{Float64}
    κ::Matrix{Float64}

    function Results(Nq::Int, Nk::Int)
        R = zeros(ComplexF64, Nq, Nk)
        R² = zeros(Float64, Nq, Nk)
        σ² = zeros(Float64, Nq, Nk)
        κ = zeros(Float64, Nq, Nk)
        return new(R, R², σ², κ)
    end
end

"""
    Preallocated(Nx, Nq, Nk)
    Preallocated(params::Parameters)

Preallocated is a struct which contains all preallocated matrices which are used in surface generation and the solution of the system of equations.
The array values of the members are mutated in place, thus each thread should have its own copy of this struct.
"""
struct Preallocated
    PMpq::Matrix{ComplexF64}
    PNpk::Matrix{ComplexF64}

    SMpq::Matrix{ComplexF64}
    SNpk::Matrix{ComplexF64}

    Fys::Vector{ComplexF64}
    sFys::Vector{ComplexF64}
    ys::Vector{Float64}

    P_res::Results
    S_res::Results

    function Preallocated(Nx::Int, Nq::Int, Nk::Int)::Preallocated
        PMpq = zeros(ComplexF64, Nq, Nq)
        PNpk = zeros(ComplexF64, Nq, Nk)
        
        SMpq = similar(PMpq)
        SNpk = similar(PNpk)

        ys = Vector{Float64}(undef, Nx)
        Fys = similar(ys, ComplexF64)
        sFys = similar(Fys)

        P_res = Results(Nq, Nk)
        S_res = Results(Nq, Nk)

        new(PMpq, PNpk, SMpq, SNpk, Fys, sFys, ys, P_res, S_res)
    end
    function Preallocated(params::Parameters)::Preallocated
        SimPrealloc(length(params.xs), length(params.qs), length(params.ks))
    end
end

"Iteratively update observable like so: ⟨A⟩ₙ = (n-1)/n ⟨A⟩ₙ₋₁ + Aₙ/n"
function observe(obs, x, n)
    return obs + (x - obs) / n
end

"Updates all observables in `res::Results` with the observation `x`"
function observe!(res::Results, x, n)
    R = res.R
    R² = res.R²
    σ² = res.σ²
    κ = res.κ

    @inbounds for I in eachindex(R)
        R[I] = observe(R[I], x[I], n)
    end

    @inbounds for I in eachindex(R²)
        R²[I] = observe(R²[I], abs2(x[I]), n)
    end

    @inbounds for I in eachindex(σ²)
        var = abs2(x[I] - R[I])
        σ²[I] = observe(σ²[I], var, n)
        κ[I] = observe(κ[I], var^2, n)
    end
end

"""
Convenience version of the above.
Updates both `pre.P_res` and `pre.S_res` with linear system solutions stored in `pre.PNpk` and `pre.SNpk`
"""
function observe!(pre::Preallocated, n::Int)
    observe!(pre.P_res, pre.PNpk, n)
    observe!(pre.S_res, pre.SNpk, n)
end


function combine(res1::Results, res2::Results, n1::Int, n2::Int)
end

"""
Combines all values in `results` weighted by the number of observations given in `counts`
"""
function combine(results::Vector{Results}, counts::Vector{Int})::Results
    out = deepcopy(results[1])
    @assert length(results) == length(counts)
    for i in eachindex
    end
end