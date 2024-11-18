

struct SolverData{ParametersType}
    params::ParametersType
    precomputed::Precomputed
    P_res::Results
    S_res::Results
    iters::Int64

    function SolverData(params::ParametersType, iters::Int64=1000) where {ParametersType}
        P_res = Results(params)
        @debug "SolverData init: Allocation for P-polarization Results completed"
        S_res = Results(params)
        @debug "SolverData init: Allocation for P-polarization Results completed"
        precomputed = Precomputed(params)
        @debug "SolverData init: Allocation for Precomputed values completed"

        return new{ParametersType}(params, precomputed, P_res, S_res, iters)
    end
end

function precompute!(data::SolverData{params})::Nothing where {params}
    precompute!(data.precomputed, data.params)
    validate(data.precomputed)
    return nothing
end
