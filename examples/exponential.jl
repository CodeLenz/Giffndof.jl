using Giffndof
include("common/base_data.jl")
include("common/process.jl")
include("common/plot.jl")


 # Run the test
 function Run_exponential(tspan  = (0.0,10.0), dt=0.001, beta_c=1E-2)

    # Process the solutions
    T,Y,U = Process( tspan, dt, Example_exponential, f_exponential!, beta_c)

    # Make the plot
    Make_plot(T,Y,U,"exponential")

 end

# Loading for Newmark
function f_exponential!(t,F::Vector{T}) where T
  
    F[1] = 0.0
    F[2] = 3*sin(4*t)
    F[3] = 0.0

end


function Example_exponential(M,C,K,U0,V0,t0)


    # Convert 3sin(4t) to exponentials
    data = Split_sin(3.0,4.0,0.0)

    # Create a dictionary. Each key corresponds to the DOF (j)
    load_data = Dict{Int64,Vector{ComplexF64}}()

    # For our example, DOF=2 and we have two exponentials
    load_data[2] = data

    # Main function -> solve the problem
    y, yh, yp = Solve_exponential(M,C,K,U0,V0,load_data,t0=t0)

    # Return the solutions for any t
    return y, yh, yp
    
 end   




