using Giffndof
using OrderedCollections
include("common/base_data.jl")
include("common/process.jl")
include("common/plot.jl")


 # Run the test
 function Run_polynomial(tspan  = (0.0,10.0), dt=0.001, beta_c=1E-6)

    # Process the solutions
    T,Y,U = Process( tspan, dt, Example_polynomial, f_polynomial!,Problem_Data,  beta_c)

    # Make the plot
    Make_plot(T,Y,U,"polynomial")

 end


# Load for Newmark
function f_polynomial!(t,F::Vector{T}) where T
  
    F[1] = 0.0
    F[2] = 10*t - t^2
    F[3] = 0.0

end



function Example_polynomial(M,C,K,U0,V0,t0)

    # Create a dictionary. Each key corresponds to the DOF (j)
    load_data = OrderedDict{Int64,Vector{Float64}}()

    # 10t - t^2 at DOF 2
    #        DOF     t2   c20  c21    c22  
    load_data[2] = [0.0 ; 0.0; 10.0; -1.0]

    # Main function -> solve the problem
    y, yh, yp = Solve_polynomial(M,C,K,U0,V0,load_data,t0=t0)

    # Return the solutions for any t
    return y, yh, yp
    
 end   


 println("You can run the example with Run_polynomial(tspan  = (0.0,10.0), dt=0.001, beta_c=1E-2)")


