
function dirac_particularResponse(times, n_times, data, response,tol_complex=1E-6)
   
       # Initializes a counter of the number of Dirac's deltas that were 
       # already turned on
       n_deltasOn = 1

       # Recover dimen
       dimen = data.dimen

       # Recover the number of terms
       nk = data.nk

       # Recover expF11
       expF11 = data.expF11_delta

       # Recover F11
       F11 = data.F11
                
       # Recover v_j
       v_j = data.v_j

       # If C_bar-F11 is complex conjugate of F11   
       if data.flag_conjugacy
            
           # Initializes a to store the response at a time instance
           response_instance = zeros(ComplexF64, dimen)

           # Iterates through the time points
           for i=1:n_times
   
               # Multiplies the vector of response at this instant by the exponential matrix
               response_instance = expF11*response_instance
   
               # Iterates through the Dirac's deltas' coefficients
               for k=n_deltasOn:nk
   
                   # Recover t_j
                   t_jk = data.t_jk[k]

                   # Recover coefficients c_j
                   c_jk = data.c_jk[k]
                   
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
   
           end
   
       # Otherwise, if they are not complex-conjugate
       else 
   
           # Initializes a to store the response at a time instance
           response_instance = zeros(ComplexF64,dimen,2)
   
           # Recover expCbF
           expCbF = data.expCF_delta

           # Recover CbF
           CbF = data.CbF

           # Iterates through the time points
           for i=1:n_times
   
               # Multiplies the vector of response at this instant by the exponential matrix
               response_instance[:,1] .= expF11*response_instance[:,1]
               response_instance[:,2] .= expCbF*response_instance[:,2]
   
               # Iterates through the Dirac's deltas' coefficients
               for k=n_deltasOn:nk
   
                   # Recover t_j
                   t_jk = data.t_jk[k]

                   # Recover coefficients c_j and v_j
                   c_jk = data.c_jk[k]
   
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

           end
   
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
   