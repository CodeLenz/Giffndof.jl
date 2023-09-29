# Defines a function to solve the particular response due to Dirac's 
# delta when Cbf and F11 are complex conjugate
function dirac_particularResponseConjugate(times, n_times, data, response,
 C2, F11, expF11, tol_complex=1E-6)
   
    # Initializes a counter of the number of Dirac's deltas that were 
    # already turned on
    n_deltasOn = 1

    # Recover dimen
    dimen = data.dimen

    # Recover the number of terms
    nk = data.nk 

    # Recover the matrix of v_j vectors
    matrix_vj = data.v_j

    # Recovers the indexes for the v_j vectors in the matrix of v_j vec-
    # tors
    indexes_vj = data.indexes_vj

    # Initializes a vector to store the response at a time instance. I-
    # nitializes it with the vector due to initial conditions C2 to in-
    # corporate the homogeneous solution in the same loop
    response_instance = C2

    # Iterates through the time points
    for i=1:n_times

        # Iterates through the Dirac's deltas' coefficients
        for k=n_deltasOn:nk

            # Recover t_j
            t_jk = data.t_jk[k]

            # Recover coefficients c_j
            c_jk = data.c_jk[k]

            # Recover the vector v_j
            v_j = matrix_vj[:,indexes_vj[k]]
            
            # If the Dirac impulse was before the current time point
            if t_jk<times[i]
                
                # As the impulse was out of any time point, its ex-
                # ponential is already live, thus, this exponential
                # term must be calculated

                # Adds the vector v_j to the response
                response_instance .+= c_jk*exp(F11*(t_jk-times[i]))*v_j

                # Updates the counter of deltas turned on
                n_deltasOn += 1

            # Otherwise, if it is exactly on the time point
            elseif t_jk==times[i]

                # Adds the vector v_j to the response
                response_instance .+= c_jk*v_j

                # Updates the counter of deltas turned on
                n_deltasOn += 1

            # The deltas must be ordered, so, if the delta is not
            # turned on, stops the search
            else

                break

            end

        end

        # Adds the instant response to the matrix of response
        @. response[:,i] += 2*real.(response_instance)

        # Multiplies the vector of response at this instant by the exponential matrix
        response_instance = expF11*response_instance
    end

    # Tests if the complex part of the response is greater than the tolerance
    complex_part = maximum(abs.(imag.(response)))

    if complex_part>tol_complex
        println("WARNING: imag(response) > tolerance",complex_part," ",tol_complex)
    end

    # Returns the response
    return nothing
   
end

# Defines a function to solve the particular response due to Dirac's 
# delta when Cbf and F11 are NOT complex conjugate
function dirac_particularResponseNotConjugate(times, n_times, data, response,
 C2, C1, F11, CbF, expF11, expCbF, tol_complex=1E-6)
   
    # Initializes a counter of the number of Dirac's deltas that were 
    # already turned on
    n_deltasOn = 1

    # Recover dimen
    dimen = data.dimen

    # Recover the number of terms
    nk = data.nk
            
    # Recover the matrix of v_j vectors
    matrix_vj = data.v_j

    # Recovers the indexes for the v_j vectors in the matrix of v_j vec-
    # tors
    indexes_vj = data.indexes_vj

    # Initializes a to store the response at a time instance. Initiali-
    # zes it with the vectors due to initial conditions, C2 and C1, to 
    # incorporate the homogeneous solution in the same loop. Adds C1 as
    # negative for, in line 164, the part relative to C1 is subtracted
    # from the part relative to C2; while in the homogeneous response,
    # this sign is omitted
    response_instance = [C2 (-C1)]

    # Iterates through the time points
    for i=1:n_times

        # Iterates through the Dirac's deltas' coefficients
        for k=n_deltasOn:nk

            # Recover t_j
            t_jk = data.t_jk[k]

            # Recover coefficients c_j and v_j
            c_jk = data.c_jk[k]

            # Recover the vector v_j
            v_j = matrix_vj[:,indexes_vj[k]]

            # Tests if the current delta was turned on already
            # If the Dirac impulse was before the current time point
            if t_jk<times[i]

                # As the impulse was out of any time point, its ex-
                # ponential is already live, thus, this exponential
                # term must be calculated

                # Adds the vector v_j to the response
                response_instance[:,1] .+= c_jk*exp(F11*(t_jk-times[i]))*v_j
                response_instance[:,2] .+= c_jk*exp(CbF*(t_jk-times[i]))*v_j

                # Updates the counter of deltas turned on
                n_deltasOn += 1

            # Otherwise, if it is exactly on the time point
            elseif t_jk==times[i]

                # Adds the vector v_j to the response
                response_instance[:,1] .+= c_jk*v_j
                response_instance[:,2] .+= c_jk*v_j

                # Updates the counter of deltas turned on
                n_deltasOn += 1

            # The deltas must be ordered, so, if the delta is not turned on, stops the search
            else

                break

            end

        end

        # Adds the instant response to the matrix of response
        @. response[:,i] += real.(response_instance[:,1]-response_instance[:,2])

        # Multiplies the vector of response at this instant by the exponential matrix
        response_instance[:,1] .= expF11*response_instance[:,1]
        response_instance[:,2] .= expCbF*response_instance[:,2]

    end

    # Tests if the complex part of the response is greater than the tolerance
    complex_part = maximum(abs.(imag.(response)))

    if complex_part>tol_complex
        println("WARNING: imag(response) > tolerance",complex_part," ",tol_complex)
    end

    # Returns the response
    return nothing
   
end
   
# Defines a function to evaluate the derivative of the particular solu-
# tion due to Dirac's delta excitation using the data in a dictionary
function dirac_derivativeParticular(t, data_dictionary, dy_p)

    return nothing

end
