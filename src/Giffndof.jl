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
    include("heaviside/commonh.jl")

    # Load solution methods
    include("exponential/solve_exponential.jl")
    include("polynomial/solve_polynomial.jl")
    include("dirac/solve_dirac.jl")
    include("heaviside/order2/solve_heaviside2.jl")
    include("heaviside/order1/solve_heaviside1.jl")
    include("heaviside/order0/solve_heaviside0.jl")
    include("heaviside_series1/solve_hs1.jl")

    # Export the methods - Solvers
    export Solve_exponential
    export Solve_polynomial
    export Solve_dirac
    export solve_heaviside2
    export solve_heaviside1
    export solve_heaviside0
    export Solve_HS1

    # Export some auxiliary functions 
    export Evaluate_coefs_c, Evaluate_gtilde

end #module
