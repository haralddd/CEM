struct Precomputed
    PMpqn::Array{ComplexF64,3}
    PNpkn::Array{ComplexF64,3}

    SMpqn::Array{ComplexF64,3}
    SNpkn::Array{ComplexF64,3}

    function Precomputed(Nq::Int, Nk::Int, Ni::Int)::Precomputed
        PMpqn = Array{ComplexF64, 3}(undef, Nq, Nq, Ni)
        PNpkn = Array{ComplexF64, 3}(undef, Nq, Nk, Ni)
        
        SMpqn = similar(PMpqn)
        SNpkn = similar(PNpkn)
        return new(PMpqn, PNpkn, SMpqn, SNpkn)
    end
    function Precomputed(params::Parameters{_S,_A,_B})::Precomputed where {_S,_A,_B}
        return Precomputed(length(params.qs), length(params.ks), params.Ni+1)
    end
    function Precomputed(params::Parameters{_S,Vacuum,Uniaxial})::Precomputed where {_S}
        return Precomputed(2*length(params.qs), length(params.ks), params.Ni+1)
    end
end

function precompute!(pre::Precomputed, params::Parameters{S,A,B})::Nothing where {S,A,B}
    M_invariant!(pre.PMpqn, params, :p)
    N_invariant!(pre.PNpkn, params, :p)

    M_invariant!(pre.SMpqn, params, :s)
    N_invariant!(pre.SNpkn, params, :s)
    return nothing
end

function _validate_single(A)
    if !all(isfinite.(A))
        open("trace.dat", "w+") do io
            write(io, A)
        end
        error("Validation of precomputed data failed, matrices not finite")
    end
end

function validate(pre::Precomputed)::Nothing
    _validate_single(pre.PMpqn)
    _validate_single(pre.PNpkn)
    _validate_single(pre.SMpqn)
    _validate_single(pre.SNpkn)
    return nothing
end