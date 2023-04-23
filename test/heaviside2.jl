#
# Test Heaviside2 solution
#
@testset "heaviside2" begin    

    # Data
    tspan  = (0.0,5.0)
    dt = 1E-4
    beta_c = 1E-3
    rtol = 1E-2

    # Loading for Newmark
    function f_heaviside2!(t,F::Vector{T}) where T
  
        F[1] = (-1+t)*Giffndof.Heaviside(t,1.0) + (1-t+t^2)*Giffndof.Heaviside(t,2.0)
        F[2] = (-t)*Giffndof.Heaviside(t,2.0) + (4+2*t-t^2)*Giffndof.Heaviside(t,3.0)
        F[3] = (-1+10*t-2*t^2)*Giffndof.Heaviside(t,4.0)

    end


    function Example_heaviside2(M,C,K,U0,V0,t0)

        # Create a dictionary. Each key corresponds to the DOF (j)
        load_data = OrderedDict{Int64,Vector{Float64}}()
    
        load_data[1] = [-1.0;1.0;0.0;1.0;  1.0;-1.0;1.0;2.0]
        load_data[2] = [0.0;-1.0;0.0;2.0;  4.0; 2.0;-1.0; 3.0]
        load_data[3] = [-1.0;10.0;-2.0;4.0]
            
        # Main function -> solve the problem
        y, yh, yp = Solve_heaviside2(M,C,K,U0,V0,load_data,t0=t0)
    
        # Return the solutions for any t
        return y, yh, yp
        
    end   
    
    # Run the test
    function Run_heaviside2(tspan, dt, beta_c, PData::Function )

        # Process the solutions
        T,Y,U = Process( tspan, dt, Example_heaviside2, f_heaviside2!, PData, beta_c)

    end

    
    function Problem_Data(beta_c=1E-6 ; U0  = [0.0; 0.0; 0.0], 
                         V0  = [0.0; 0.0; 0.0])

        # Mass matrix
        M = [2.0 0.0 0.0 ;
            0.0 2.0 0.0 ;
            0.0 0.0 1.0 ]

        # Stiffness matrix
        K = [6.0 -4.0  0.0 ;
            -4.0  6.0 -2.0 ;
            0.0 -2.0  6.0]*1E2

        # Damping matrix
        C = beta_c*K

        # at time t0
        t0 = 0.0

        # Return the data
        return M,C,K,U0,V0,t0

    end


    # Call the solution procedures - ours and Newmark
    T,Y,U = Run_heaviside2(tspan, dt, beta_c, Problem_Data)

    # Test if T is ok
    @test isapprox(T,collect(tspan[1]:dt:tspan[2]),rtol=rtol)

    # Test if the DOFs are OK by comparing Y to U
    for i=1:3
        @test isapprox(Y[i,:],U[:,i],rtol=rtol)
    end

    # The same for non homogeneous IC 

    # Call the solution procedures - ours and Newmark
    PD(beta_c) = Problem_Data(beta_c,U0=[0.1;-0.1;0.3])
    T,Y,U = Run_heaviside2(tspan, dt, beta_c, PD)

    # Test if the DOFs are OK by comparing Y to U
    for i=1:3
        @test isapprox(Y[i,:],U[:,i],rtol=rtol)
    end

end # testset






