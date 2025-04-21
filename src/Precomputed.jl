"""
Precomputed struct for Hybrid solver
    Contains the elements required to calculate T from R, which supplements the base precomputed struct
    in the Hybrid solver.
"""
struct HybridPrecomputed
    PAqn::Matrix{ComplexF64}
    PB1kn::Matrix{ComplexF64}
    PB2qn::Matrix{ComplexF64}

    SAqn::Matrix{ComplexF64}
    SB1kn::Matrix{ComplexF64}
    SB2qn::Matrix{ComplexF64}

    function HybridPrecomputed(Nq::Int, Nk::Int, Ni::Int)
        PAqn = Matrix{ComplexF64}(undef, Nq, Ni)
        PB1kn = Matrix{ComplexF64}(undef, Nk, Ni)
        PB2qn = Matrix{ComplexF64}(undef, Nq, Ni)
        SAqn = Matrix{ComplexF64}(undef, Nq, Ni)
        SB1kn = Matrix{ComplexF64}(undef, Nk, Ni)
        SB2qn = Matrix{ComplexF64}(undef, Nq, Ni)
        return new(PAqn, PB1kn, PB2qn, SAqn, SB1kn, SB2qn)
    end
    function HybridPrecomputed(data::SolverData)
        Nq = length(data.params.qs)
        Nk = length(data.params.ks)
        Ni = data.params.Ni+1
        return HybridPrecomputed(Nq, Nk, Ni)
    end
end

struct Precomputed
    PMpqn::Array{ComplexF64,3}
    PNpkn::Array{ComplexF64,3}

    SMpqn::Array{ComplexF64,3}
    SNpkn::Array{ComplexF64,3}

    hybrid::Union{HybridPrecomputed, Nothing}

    function Precomputed(Nq::Int, Nk::Int, Ni::Int, hybrid::Union{HybridPrecomputed, Nothing})::Precomputed
        PMpqn = Array{ComplexF64, 3}(undef, Nq, Nq, Ni)
        PNpkn = Array{ComplexF64, 3}(undef, Nq, Nk, Ni)
        
        SMpqn = similar(PMpqn)
        SNpkn = similar(PNpkn)
        return new(PMpqn, PNpkn, SMpqn, SNpkn, hybrid)
    end
    function Precomputed(data::SolverData)::Precomputed
        if data.solver_type == :full
            return Precomputed(2*length(data.params.qs), length(data.params.ks), data.params.Ni+1, nothing)
        elseif data.solver_type == :hybrid
            return Precomputed(length(data.params.qs), length(data.params.ks), data.params.Ni+1, HybridPrecomputed(data))
        else
            return Precomputed(length(data.params.qs), length(data.params.ks), data.params.Ni+1, nothing)
        end
    end
end

function validate(pre::HybridPrecomputed)::Nothing
    @assert all(isfinite.(pre.PAqn))
    @assert all(isfinite.(pre.PB1kn))
    @assert all(isfinite.(pre.PB2qn))

    @assert all(isfinite.(pre.SAqn))
    @assert all(isfinite.(pre.SB1kn))
    @assert all(isfinite.(pre.SB2qn))
    return nothing
end

function validate(pre::Precomputed)::Nothing
    @assert all(isfinite.(pre.PMpqn))
    @assert all(isfinite.(pre.PNpkn))
    @assert all(isfinite.(pre.SMpqn))
    @assert all(isfinite.(pre.SNpkn))
    if pre.hybrid !== nothing
        validate(pre.hybrid)
    end
    return nothing
end

function precompute!(hybrid::HybridPrecomputed, data::SolverData)::Nothing
    A_invariant_hybrid!(hybrid.PAqn, data.params, :p)
    B1_invariant_hybrid!(hybrid.PB1kn, data.params)
    B2_invariant_hybrid!(hybrid.PB2qn, data.params)
    A_invariant_hybrid!(hybrid.SAqn, data.params, :s)
    B1_invariant_hybrid!(hybrid.SB1kn, data.params)
    B2_invariant_hybrid!(hybrid.SB2qn, data.params)
    return nothing
end

function precompute!(pre::Precomputed, data::SolverData)::Nothing
    if data.solver_type == :full
        M_invariant_full!(pre.PMpqn, data.params, :p)
        N_invariant_full!(pre.PNpkn, data.params, :p)
        M_invariant_full!(pre.SMpqn, data.params, :s)
        N_invariant_full!(pre.SNpkn, data.params, :s)
    else
        M_invariant_reduced!(pre.PMpqn, data.params, :p)
        N_invariant_reduced!(pre.PNpkn, data.params, :p)
        M_invariant_reduced!(pre.SMpqn, data.params, :s)
        N_invariant_reduced!(pre.SNpkn, data.params, :s)
        if pre.hybrid !== nothing
            precompute!(pre.hybrid, data)
        end
    end
    validate(pre)
    return nothing
end