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

    function Preallocated(Nx::Int, Nq::Int, Nk::Int)::Preallocated
        PMpq = zeros(ComplexF64, Nq, Nq)
        PNpk = zeros(ComplexF64, Nq, Nk)
        
        SMpq = zeros(ComplexF64, Nq, Nq)
        SNpk = zeros(ComplexF64, Nq, Nk)

        ys = Vector{Float64}(undef, Nx)
        Fys = similar(ys, ComplexF64)
        sFys = similar(Fys)

        new(PMpq, PNpk, SMpq, SNpk, Fys, sFys, ys)
    end
    function Preallocated(params::Parameters)::Preallocated
        Preallocated(length(params.xs), length(params.qs), length(params.ks))
    end
    function Preallocated(params::Parameters{_S,Vacuum,Uniaxial}) where {_S}
        Preallocated(length(params.xs), 2*length(params.qs), length(params.ks))
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