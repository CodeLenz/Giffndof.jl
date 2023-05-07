#
# Test HS1 solution
#
@testset "HS1" begin    

    # Data
    tspan  = (0.0,5.0)
    dt = 1E-4
    beta_c = 1E-3
    rtol = 1E-2

    # Loading
    function g(t)
        -cos(0.5*t) + sin(t)  + cos(1.5*t - 1.5) -2*sin(2*t) + 2*sin(10*t)  
    end

    # Load for Newmark
    function f_hs1!(t,F::Vector{T}) where T
    
        F[1] = 0
        F[2] = g(t)
        F[3] = 0

    end


    function Example_hs1(M,C,K,U0,V0,t0)

        # Loading
        load_data = OrderedDict{Int64,Function}()

        # Pass function g(t) to the dictionary
        load_data[2] = g

        # Create Ts by using tspan and dt
        Ts = tspan[1]:dt:tspan[2]

        #  Main function -> solve the problem
        y, yh, yp = Solve_HS1(M,C,K,U0,V0,Ts,load_data,t0=t0)
    
        # Return the solutions for any t
        return y, yh, yp
        
    end   
    

    # Run the test
    function Run_hs1(tspan, dt, beta_c, PData::Function )

        # Process the solutions
        T,Y,U = Process( tspan, dt, Example_hs1, f_hs1!, PData, beta_c)

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
    T,Y,U = Run_hs1(tspan, dt, beta_c, Problem_Data)

    # Test if T is ok
    @test isapprox(T,collect(tspan[1]:dt:tspan[2]),rtol=rtol)

    # Test if the DOFs are OK by comparing Y to U
    for i=1:3
        @test isapprox(Y[i,:],U[:,i],rtol=rtol)
    end

    # The same for non homogeneous IC 

    # Call the solution procedures - ours and Newmark
    PD(beta_c) = Problem_Data(beta_c,U0=[0.1;-0.1;0.3])
    T,Y,U = Run_hs1(tspan, dt, beta_c, PD)

    # Test if the DOFs are OK by comparing Y to U
    for i=1:3
        @test isapprox(Y[i,:],U[:,i],rtol=rtol)
    end

end # testset


