# Routine to store the structs for the particular response functions

# Exponential
struct exponential_dataStruct{T}

    # Defines a counter of exponential terms
    nk::Int64

    # Defines the vectors of coefficients c_jk, beta_jk and phi_jk of the exponential terms
    c_jk::Vector{T}

    beta_jk::Vector{T}

    phi_jk::Vector{T}

    # Defines a matrix to store the k_jk vectors at each column
    matrix_kjk::Matrix{T}

end

# Polynomial
struct polynomial_dataStruct{T,C}

    # Defines a counter of polynomial terms
    nk::Int64

    # Defines the vectors of coefficients c_jk, time shifts t_j and powers k of the polynomial terms
    c_j::Vector{T}

    t_j::Vector{T}

    degrees::Vector{Int64}

    # Defines a matrix to store the v_k_l_p ([M*(F11^(l-p))*(K_bar^p)]\e_j) vectors at each column
    matrix_vklp::Matrix{C}

end

# Dirac's delta when F11 and C_bar-F11 are complex-conjugate
struct dirac_dataStructConjugate{T,C}

    # Defines the problem's dimensionality
    dimen::Int64

    # Defines a counter of delta terms
    nk::Int64

    # Defines the vectors of coefficients c_jk and time shifts t_jk
    c_jk::Vector{T}

    t_jk::Vector{T}

    # Defines a vector to store the vector ((C-(2*M*F11))\input_vector)
    v_j::Vector{C}

    # Defines a flag to signal whether F11 and C_bar-F11 are complex-conjugates or not
    flag_conjugacy::Bool

    # Defines a matrix to F11
    F11::Matrix{C}

    # Defines a matrix to store the exponentials
    expF11_delta::Matrix{C}

end

# Dirac's delta when F11 and C_bar-F11 are not complex-conjugate
struct dirac_dataStructNonConjugate{T,C}

    # Defines the problem's dimensionality
    dimen::Int64

    # Defines a counter of delta terms
    nk::Int64

    # Defines the vectors of coefficients c_jk and time shifts t_jk
    c_jk::Vector{T}

    t_jk::Vector{T}

    # Defines a vector to store the vector ((C-(2*M*F11))\input_vector)
    v_j::Vector{C}

    # Defines a flag to signal whether F11 and C_bar-F11 are complex-conjugates or not
    flag_conjugacy::Bool

    # Defines a matrix to F11
    F11::Matrix{C}

    # Defines a matrix to C_bar-F11
    CbF::Matrix{C}

    # Defines matrix exponentials
    expF11_delta::Matrix{C}
    expCF_delta::Matrix{C}

end