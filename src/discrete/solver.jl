#
#
#
# Defines a function to evaluate the integration constants
#
#
"""
   Solve the system of second order ODEs 

   ``M A(t) + C V(t) + K U(t) = F(t)``

   with 

   U(0.0) = U0

   V(0.0) = V0

   at DISCRETE times.

   Loading is informed by using a dictionary  ``load_data::Dict{Int64, Matrix{Union{String,Float64}}}``

   The output is an n x nt matrix where the rows are the DOFs and the columns the values
   for each discrete time.

"""
function Solve_discrete(M::AbstractMatrix{TF},C::AbstractMatrix{TF},K::AbstractMatrix{TF},
                        times::Ts,load_data::Dict{Int64, Matrix{Union{String,Float64}}},
                        U0::Vector{TF},V0::Vector{TF}; tol=1E-6) where {Ts,TF}

    # Number of time points
    n_times = length(times)

    # Test if we have at least two time points
    n_times>2 || throw("Solve::we need at least two time steps")

    # Calculates the time delta
    dt = times[2]-times[1]

    # Test if dt is positive 
    dt>=0 || throw("Solve::time step must be >=0")

    # Initial time
    t0 = times[1]

    # Asserts that all the matrices have the same dimension
    size(M)==size(K)==size(C) || throw("Solve::input matrices must have the same dimensions") 

    # Assert that matrices are square
    size(M,1)==size(M,2) || throw("Solve::input matrices must be square") 

    # Obtains the dimension of the problem (number of DOFs)
    dimen = size(M,1)

    # Create the vector of excited dofs
    excited_dofs = collect(keys(load_data))
    sort!(excited_dofs)
    unique!(excited_dofs)

    # Number of excited DOFs
    n_excitedDOF = length(excited_dofs)

    # Test if we have at least one loaded DOF
    n_excitedDOF>=1 || through("Solve::we need at least a loaded DOF")

    # Asssert that excited_dofs are within the bounds 
    maximum(excited_dofs)<=dimen || throw("Solve::excited_dofs must be within the number of dofs of the problem")
    minimum(excited_dofs)>=1     || throw("Solve::excited_dofs must be within the number of dofs of the problem")
    
    # Factorizes the mass matrix
    choleskyM = cholesky(Symmetric(M))

    # Calculates C_bar and K_bar
    C_bar = Array(choleskyM\C)
    K_bar = Array(choleskyM\K)

    # Calculates F211
    F11 = 0.5*(C_bar.+sqrt((C_bar^2).-(4*K_bar)))

    # Calculates the exponential of F11 multiplied by a delta
    expF11_delta = exp(-F11*dt)

    # Calculates CbF 
    CbF = C_bar .- F11

    # Initializes the matrix of the response
    response = zeros(ComplexF64, dimen, n_times)

    # Initializes a vector of the derivative of the particular response at time t0
    dy_p = zeros(ComplexF64, dimen)

    # Creates the input vector
    input_vector = zeros(dimen)

    # Initializes the data structure for Dirac excitation and checks if
    # F11 is complex
    
    complex_partF11 = norm(imag(F11))
    
    if complex_partF11>0

        data_dirac = dirac_dataStruct(dimen, 0, Float64[], Float64[], 
        zeros(ComplexF64, dimen, 1), Int64[])

    else

        # Initializes the data structure for Dirac excitation
        data_dirac = dirac_dataStruct(dimen, 0, Float64[], Float64[], 
        zeros(Float64, dimen, 1), Int64[])

    end

    # Iterates through the excited DOFs
    for i=1:n_excitedDOF

        # Set DOF to 1.0
        input_vector[excited_dofs[i]] = 1.0

        # Build the dictionary of data for the particular solution
        data_exponential, data_polynomial, data_dirac, flag_exp, flag_pol = Process(
         dimen,M,C,K,load_data[excited_dofs[i]], F11, C_bar, K_bar,
         input_vector, tol, data_dirac)

        # If the exponential particular solution is required
        if flag_exp 

            # Updates the matrix of response
            exponential_particularResponse(times,n_times, data_exponential, response, tol)

            # Updates the vector of the derivative of the response at time t0
            exponential_derivativeParticular(t0,data_exponential,dy_p)

        end

        # If the polynomial particular solution is required
        if flag_pol 

            # Updates the matrix of response
            polynomial_particularResponse(times,n_times,data_polynomial,response, tol)

            # Updates the vector of the derivative of the response at time t0
            polynomial_derivativeParticular(t0,data_polynomial,dy_p)

        end

        # UnSet the DOF
        input_vector[excited_dofs[i]] = 0.0

    end

    # Calculates integrating constants for homogeneous initial condi-
    # tions at t0=0
    C1, C2 = calculate_integrationConstants(C_bar, F11, response[:,1], dy_p, U0, V0)

    # Adds the response due to Dirac's delta and the homogeneous solu-
    # tion. If there is no Dirac's delta excitation, adds the homoge-
    # neous response only

    # Check the conjugacy between CbF and F11
    norm_conjugacy = norm((CbF).-conj(F11))
   
    if norm_conjugacy<tol

        dirac_particularResponseConjugate(times, n_times, data_dirac,
         response, C2, F11, expF11_delta, tol)

    else

        # Calculates the exponential of Cb-F11
        expCF_delta = exp(-dt*CbF)

        dirac_particularResponseNotConjugate(times, n_times, data_dirac,
         response, C2, C1, F11, CbF, expF11_delta, expCF_delta, tol)

    end

    # Returns the response
    return response

end