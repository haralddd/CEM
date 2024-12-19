function chunk(i)
    return "Hello from $i"
end

tasks = fetch.([Threads.@spawn chunk(i) for i in 1:Threads.nthreads()])