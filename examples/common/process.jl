function Process(tspan::Tuple{Float64,Float64}, dt::Float64, Example::Function, 
                 f!::Function, FData::Function,
                 beta_c = 1E-2)

    # Generate the discrete times
    tt = tspan[1]:dt:tspan[2]

    # Load Problem data
    M,C,K,U0,V0,t0 = FData(beta_c)

    # Solve the problem
    y, _ =  Example(M,C,K,U0,V0,t0)

    # Reshape the outputs in a discrete grid
    ndofs = size(y(0.0),1)
    Y = reshape([real(y(t))[k] for k=1:ndofs for t in tt],length(tt),ndofs)

    # Solve by using Newmark - Beta method
    U,T = Solve_newmark(M, C, K, f!, tspan, dt, U0=U0, V0=V0)

    # Return the solutions
    return T, Y, U    
    

end
