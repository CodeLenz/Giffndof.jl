using Giffndof
using OrderedCollections

include("common/base_data.jl")
include("common/process.jl")
include("common/plot.jl")


 # Run the test
 function Run_heaviside2(tspan  = (0.0,10.0), dt=0.001, beta_c=1E-2)

    # Process the solutions
    T,Y,U = Process( tspan, dt, Example_heaviside2, f_heaviside2!, beta_c)

    # Make the plot
    Make_plot(T,Y,U,"heaviside2")

 end


# Loading for Newmark
function f_heaviside2!(t,F::Vector{T}) where T
  
    val  = (-30+40*t-10*t^2)*Giffndof.Heaviside(t,1.0) 
    val += ( 30-40*t+10*t^2)*Giffndof.Heaviside(t,3.0) 

    F[1] = 0.0
    F[2] = val
    F[3] = 0.0

end


function Example_heaviside2(M,C,K,U0,V0,t0)

    #
    # Loading (-30 + 40*t - 10*t^2) H(t-1) + (30 - 40*t + 10*t^2)H(t-3)
    #
    load_data = OrderedDict{Int64,Vector{Float64}}()

    #   c_j00 c_j01 c_j02  t_jk .... c_j(nk)0 c_j(nk)1 c_j(nk)2 t_j(nk)
    load_data[2] = [-30.0; 40.0; -10.0; 1.0; 30.0; -40.0; 10.0; 3.0]
    
    # Main function -> solve the problem
    y, yh, yp = Solve_heaviside2(M,C,K,U0,V0,load_data,t0=t0)

    # Return the solutions for any t
    return y, yh, yp
    
 end   


 println("You can run the example with Run_heaviside2(tspan  = (0.0,10.0), dt=0.001, beta_c=1E-2)")
