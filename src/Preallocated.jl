"""
    Results(Nq, Nk)
    Results(params::Parameters)

Results is a struct which contains all preallocated matrices which are used in surface generation and the solution of the system of equations.
The array values of the members are mutated in place, thus each thread should have its own copy of this struct.
"""
struct Results
    A::Matrix{ComplexF64}
    A²::Matrix{Float64}
    σ²::Matrix{Float64}
    κ::Matrix{Float64}

    function Results(Nq::Int, Nk::Int)
        A = zeros(ComplexF64, Nq, Nk)
        A² = zeros(Float64, Nq, Nk)
        σ² = zeros(Float64, Nq, Nk)
        κ = zeros(Float64, Nq, Nk)
        return new(A, A², σ², κ)
    end
    function Results(params::Parameters)
        return Results(length(params.qs), length(params.ks))
    end
end

"""
    Struct for the purpose of serialization.
    Containing the results and parameters of a solver run

"""
struct SolverData{ParametersType}
    params::ParametersType
    Rp::Results
    Rs::Results
    Tp::Results
    Ts::Results
    iters::Int64
    solver_type::Symbol

    function SolverData(params::ParametersType, iters::Int64=1000, solver_type::Symbol=:full) where {ParametersType}
        Rp = Results(params)
        Rs = Results(params)
        Tp = Results(params)
        Ts = Results(params)
        @assert solver_type in [:full, :reduced, :hybrid]
        return new{ParametersType}(params, Rp, Rs, Tp, Ts, iters, solver_type)
    end
end

"""
    HybridPreallocated(Nq, Nk)
    HybridPreallocated(params::Parameters)

HybridPreallocated is a struct which contains solution vectors for the hybrid solver.
"""
struct HybridPreallocated
    PApq::Matrix{ComplexF64}
    Pbpk::Matrix{ComplexF64}

    SApq::Matrix{ComplexF64}
    Sbpk::Matrix{ComplexF64}

    function HybridPreallocated(Nq::Int, Nk::Int)
        PApq = Matrix{ComplexF64}(undef, Nq, Nq)
        Pbpk = Matrix{ComplexF64}(undef, Nq, Nk)
        SApq = Matrix{ComplexF64}(undef, Nq, Nq)
        Sbpk = Matrix{ComplexF64}(undef, Nq, Nk)
        return new(PApq, Pbpk, SApq, Sbpk)
    end
    function HybridPreallocated(data::SolverData)
        return HybridPreallocated(length(data.params.qs), length(data.params.ks))
    end
end

"""
    Preallocated(Nx, Nq, Nk)
    Preallocated(params::Parameters)

Preallocated is a struct which contains all preallocated matrices which are used in surface generation and the solution of the system of equations.
The array values of the members are mutated in place, thus each thread should have its own copy of this struct.
For the full scattering integral solver, both R and T are stored in the data members
for more efficient solving of the linear system.
"""
struct Preallocated
    PMpq::Matrix{ComplexF64}
    PNpk::Matrix{ComplexF64}

    SMpq::Matrix{ComplexF64}
    SNpk::Matrix{ComplexF64}

    Fys::Vector{ComplexF64}
    sFys::Vector{ComplexF64}
    ys::Vector{Float64}

    hybrid::Union{HybridPreallocated, Nothing}

    function Preallocated(Nx::Int, Nq::Int, Nk::Int, hybrid::Union{HybridPreallocated, Nothing})::Preallocated
        PMpq = zeros(ComplexF64, Nq, Nq)
        PNpk = zeros(ComplexF64, Nq, Nk)
        
        SMpq = zeros(ComplexF64, Nq, Nq)
        SNpk = zeros(ComplexF64, Nq, Nk)

        ys = Vector{Float64}(undef, Nx)
        Fys = similar(ys, ComplexF64)
        sFys = similar(Fys)

        new(PMpq, PNpk, SMpq, SNpk, Fys, sFys, ys, hybrid)
    end
    function Preallocated(data::SolverData)::Preallocated
        Nx = length(data.params.xs)
        Nq = length(data.params.qs)
        Nk = length(data.params.ks)
        if data.solver_type == :full
            return Preallocated(Nx, 2*Nq, Nk, nothing)
        elseif data.solver_type == :hybrid
            return Preallocated(Nx, Nq, Nk, HybridPreallocated(data))
        else
            return Preallocated(Nx, Nq, Nk, nothing)
        end
    end
end




"Iteratively update observable like so: ⟨A⟩ₙ = (n-1)/n ⟨A⟩ₙ₋₁ + Aₙ/n"
function observe(obs, x, n)
    return obs + (x - obs) / n
end

"Updates all observables in `res::Results` with the observation `x`"
function observe!(res::Results, x, n)
    A = res.A
    A² = res.A²
    σ² = res.σ²
    κ = res.κ

    for I in eachindex(A)
        A[I] = observe(A[I], x[I], n)
    end

    for I in eachindex(A²)
        A²[I] = observe(A²[I], abs2(x[I]), n)
    end

    for I in eachindex(σ²)
        var = abs2(x[I] - A[I])
        σ²[I] = observe(σ²[I], var, n)
        κ[I] = observe(κ[I], var^2, n)
    end
end

"""
    Updates all `Results` observables in `data::SolverData` with the observations from `pre::Preallocated`.
    
    # Arguments:
    - `data`: [`SolverData`](@ref) containing the results and parameters
    - `pre`: [`Preallocated`](@ref) containing preallocated values
    - `n`: Number of iterations
"""
function observe!(data::SolverData, pre::Preallocated, n::Int)
    solver_type = data.solver_type
    if solver_type == :full
        half = size(pre.PNpk, 1) ÷ 2
        Rp = @view pre.PNpk[1:half, :]
        Rs = @view pre.SNpk[1:half, :]
        Tp = @view pre.PNpk[half+1:end, :]
        Ts = @view pre.SNpk[half+1:end, :]
        observe!(data.Rp, Rp, n)
        observe!(data.Rs, Rs, n)
        observe!(data.Tp, Tp, n)
        observe!(data.Ts, Ts, n)
    elseif solver_type == :reduced
        observe!(data.Rp, pre.PNpk, n)
        observe!(data.Rs, pre.SNpk, n)
    elseif solver_type == :hybrid
        observe!(data.Rp, pre.PNpk, n)
        observe!(data.Rs, pre.SNpk, n)
        observe!(data.Tp, pre.hybrid.Pbpk, n)
        observe!(data.Ts, pre.hybrid.Sbpk, n)
    else
        error("Unknown solver type: $(solver_type)")
    end
end


"""
Adds the results `obs` to `out` with the weight `N/divisor`. 
Used to combine several means with different number of iterations, typically after multithreading
"""
function combine!(out::Results, obs::Results, N::Int, divisor::Int)
    A = out.A
    A² = out.A²
    σ² = out.σ²
    κ = out.κ

    for I in eachindex(A)
        A[I] += obs.A[I] * N / divisor
    end

    for I in eachindex(A)
        A²[I] += obs.A²[I] * N / divisor
    end

    for I in eachindex(A)
        σ²[I] += obs.σ²[I] * N / divisor
    end

    for I in eachindex(A)
        κ[I] += obs.κ[I] * N / divisor
    end
end