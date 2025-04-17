
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

    function SolverData(params::ParametersType, iters::Int64=1000) where {ParametersType}
        Rp = Results(params)
        Rs = Results(params)
        Tp = Results(params)
        Ts = Results(params)
        return new{ParametersType}(params, Rp, Rs, Tp, Ts, iters)
    end
end


function observe!(data::SolverData, pre::Preallocated, n::Int, type::Symbol=:full)
    if type == :full
        half = size(pre.PNpk, 2) ÷ 2
        Rp = @view pre.PNpk[1:half, :]
        Rs = @view pre.SNpk[1:half, :]
        Tp = @view pre.PNpk[half+1:end, :]
        Ts = @view pre.SNpk[half+1:end, :]
        observe!(data.Rp, Rp, n)
        observe!(data.Rs, Rs, n)
        observe!(data.Tp, Tp, n)
        observe!(data.Ts, Ts, n)
    elseif type == :rre
        observe!(data.Rp, pre.PNpk, n)
        observe!(data.Rs, pre.SNpk, n)
        calculate_T!(data, pre, n)
    end
end