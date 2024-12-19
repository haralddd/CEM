include("testconfig.jl")

using Plots
function test_single_solve_uniaxial()
    data = config_default_uniaxial()
    @info "precompute"
    @time precompute!(data.precomputed, data.params)
    try
        validate(data.precomputed)
    catch e
        open("trace.dat") do io
            A = similar(data.precomputed.PMpqn)
            read!(io, A)
            hms = [heatmap(.!isfinite.(A[:,:,n]), cbar=false) for n in axes(A, 3)]
            for n in axes(A, 3)
                if !all(isfinite.(A[:,:,n]))
                    println("$n not finite")
                    idxs = findall(.!isfinite.(A[:,:,n]))
                    println(idxs)
                else
                    println("$n finite")
                end
            end
            cols = 3
            rows = ceil(Int, size(A, 3) / cols)
            plot(hms..., layout=(rows,cols), legend=false, size=(cols*300,rows*300)) |> display
        end
        error("validation failed with $e")
    end

    alloc = Preallocated(data.params)
    @info "generate surface"
    @time generate_surface!(alloc, data.params)
    @info "solve single"
    @time solve_single!(alloc, data)
    plot(abs2.(alloc.SNpk[:, 1]))
end

test_single_solve_uniaxial()