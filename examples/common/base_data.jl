
function Problem_Data(beta_c=1E-6)

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

    # Initial Conditions
    U0  = [0.0; 0.0; 0.0]
    V0  = [0.0; 0.0; 0.0]

    # at time t0
    t0 = 0.0

    # Return the data
    return M,C,K,U0,V0,t0

end
