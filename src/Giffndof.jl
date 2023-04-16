module Giffndof

    # If possible, set optimization to 3
    if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
        @eval Base.Experimental.@optlevel 3
    end

    # load packages
    using LinearAlgebra, OrderedCollections, QuadGK

    # load common functions
    include("common/homo.jl")
    include("common/step.jl")

    # Load solution methods
    include("exponential/solve_exponential.jl")
    include("polynomial/solve_polynomial.jl")
    include("dirac/solve_dirac.jl")
    include("heaviside/solve_heaviside.jl")
    include("heaviside_series1/solve_hs1.jl")

    # Export the methods - Solvers
    export Solve_exponential
    export Solve_polynomial
    export Solve_dirac
    export Solve_heaviside
    export Solve_HS1

    # Export some auxiliary functions
    export Evaluate_coefs_c, Create_gtilde

end #module
