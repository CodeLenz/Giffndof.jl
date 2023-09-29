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

# Dirac's delta
mutable struct dirac_dataStruct{T,C}

    # Defines the problem's dimensionality
    dimen::Int64

    # Defines a counter of delta terms
    nk::Int64

    # Defines the vectors of coefficients c_jk and time shifts t_jk
    c_jk::Vector{T}

    t_jk::Vector{T}

    # Defines a vector to store the vector ((C-(2*M*F11))\input_vector)
    v_j::Matrix{C}

    # Defines a vector of indexes for the vectors v_j in the matrix of 
    # v_j vector

    indexes_vj::Vector{Int64}

end