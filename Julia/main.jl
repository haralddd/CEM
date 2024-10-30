include("cli.jl")
@debug PROGRAM_FILE
@debug @__FILE__
if (abspath(PROGRAM_FILE) == @__FILE__) 
    cli_main()
elseif isinteractive()
    # Do interactive stuff here
else
    error("Unknown environment")
    exit(1)
end
