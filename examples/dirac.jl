using Giffndof
using OrderedCollections
include("common/base_data.jl")
include("common/process.jl")
include("common/plot.jl")

function Example_dirac(M,C,K,U0,V0,t0)

    load_data = OrderedDict{Int64,Vector{Float64}}()

    #     g_2(t) = delta(t-1) - delta(t-5)
    #
    #               c_20  t_20  c_21  t_21
    load_data[2] = [1.0 ; 1.0; -1.0 ; 5.0]

    # Main function -> solve the problem
    y, yh, yp = Solve_dirac(M,C,K,U0,V0,load_data,t0=t0)

    # Return the solutions for any t
    return y, yh, yp
    
 end   


# Emulate dirac
function Impulse(t,t0,eps,A)

    a = 1/(2*eps)
    val = 0.0

    if t0-eps <= t <= t0+eps
        val = A*a*(1+cos(pi*(t-t0)/eps))
    end

    return val

end

# Load for Newmark
function f_dirac!(t,F::Vector{T}) where T
  
    val1 = Impulse(t,1.0,0.01,1.0)
    val5 = Impulse(t,5.0,0.01,-1.0)

    @show t, val1, val5

    F[1] = 0.0
    F[2] = val1 + val5
    F[3] = 0.0

end

 # Run the test
 function Run_dirac(tspan=(0.0,10.0), dt=0.001, beta_c=1E-2)

    # Process the solutions
    T,Y,U = Process(tspan, dt, Example_dirac, f_dirac!, beta_c)

    # Make the plot
    Make_plot(T,Y,U,"dirac")

 end
