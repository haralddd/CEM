# using LoopVectorization
using BenchmarkTools
using Statistics
using Random

observe1(obs, x, n) = obs*(n-1)/n + x/n
observe2(obs, x, n) = (obs*(n-1) + x)/n
observe3(obs, x, n) = obs + (x-obs)/n
rand_val(rng) = randn(rng) * 1.0e8
custom_rng(seed) = Random.seed!(Random.default_rng(), seed); Random.default_rng()

function rng_seq_stable(N=100, seed=1234)
    # Sequential stable mean
    rng = custom_rng(seed)
    _ret = zeros(Float64, N)
    for i in 1:N
        myrand = rand_val(rng)
        _ret[i] = myrand
    end
    return mean(_ret)
end

function rng_seq_acc(N=100, seed=1234)
    # Sequential accumulator mean
    rng = custom_rng(seed)
    ret = 0.0
    for i in 1:N
        myrand = rand_val(rng)
        ret += myrand
    end
    ret /= N
    return ret
end

function rng_seq_obs1(N=100, seed=1234)
    # Sequential observe1
    rng = custom_rng(seed)
    ret = 0.0
    for i in 1:N
        myrand = rand_val(rng)
        ret = observe1(ret, myrand, i)
    end
    return ret
end

function rng_seq_obs2(N=100, seed=1234)
    # Sequential observe2
    rng = custom_rng(seed)
    ret = 0.0
    for i in 1:N
        myrand = rand_val(rng)
        ret = observe2(ret, myrand, i)
    end
    return ret
end

function rng_seq_obs3(N=100, seed=1234)
    # Sequential observe3
    rng = custom_rng(seed)
    ret = 0.0
    for i in 1:N
        myrand = rand_val(rng)
        ret = observe3(ret, myrand, i)
    end
    return ret
end

function rng_para_direct(N=100, seed=1234)
    # Parallel direct mean
    rng = custom_rng(seed)
    ret = zeros(Float64, N)
    Threads.@threads for i in 1:N
        myrand = rand_val(rng)
        ret[i] = myrand
    end
    return mean(ret)
end

function rng_para_acc(N=100, seed=1234)
    # Parallel accumulator mean
    rng = custom_rng(seed)
    ret = Threads.Atomic{Float64}(0.0)
    Threads.@threads for i in 1:N
        myrand = rand_val(rng)
        Threads.atomic_add!(ret, myrand)
    end
    ret.value /= N
    return ret.value
end

function rng_para_direct_local(N=100, seed=1234)
    # Parallel accumulator mean
    ret = zeros(Float64, N)
    for i in 1:N
        rng = custom_rng(seed)
        myrand = rand_val(rng)
        ret[i] = myrand
    end
    ret.value /= N
    return ret.value
end

function test_rng_sequence()
    N = 1000
    M = 1000

    names = ["rng_seq_acc", "rng_seq_obs1", "rng_seq_obs2", "rng_seq_obs3",
        "rng_para_direct", "rng_para_acc"]
    funcs = [rng_seq_acc, rng_seq_obs1, rng_seq_obs2, rng_seq_obs3,
        rng_para_direct, rng_para_acc]
    diffs = zeros(length(funcs), M)

    for fi in eachindex(funcs)
        @info names[fi]
        @btime $(funcs[fi])($N, 1234)
    end

    @info "Mean errors:"
    @time for m in 1:M
        seed = rand(0:typemax(Int))

        res0 = rng_seq_stable(N, seed)

        for fi in eachindex(funcs)
            res = funcs[fi](N, seed)
            diffs[fi, m] = abs(res - res0)
        end
    end

    mean_errs = mean(diffs, dims=2)
    
    println(["$(names[i]):  $(mean_errs[i])\n" for i in eachindex(names)]...)

    nothing
end

test_rng_sequence()
