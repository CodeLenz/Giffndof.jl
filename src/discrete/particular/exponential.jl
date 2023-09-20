# Defines a function to evaluate the particular solution due to exponen-
# tial excitation using the data in a dictionary
function exponential_particularResponse(times, n_times,data,response,tol_complex=1E-6)
   
       # Takes some variables that are to be used
       c_jk    = data.c_jk
       beta_jk = data.beta_jk
       phi_jk  = data.phi_jk
   
       # Iterates through the k terms in the exponential excitation
       for k=1:data.nk
   
           # Access the vector k_jk
           k_jk = @view data.matrix_kjk[:,k]
   
           # Iterates through the time points
           for m=1:n_times   
               @. response[:,m] += c_jk[k]*exp((beta_jk[k]*times[m])+phi_jk[k])*k_jk
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
   # tion due to exponential excitation using the data in a dictionary
   function exponential_derivativeParticular(t, data, dy_p)
   
       
       # Takes some variables that are to be used
       c_jk    = data.c_jk
       beta_jk = data.beta_jk
       phi_jk  = data.phi_jk
   
       # Iterates through the k terms in the exponential excitation
       for k=1:data.nk
   
           # Access the vector k_jk
           k_jk = @view data.matrix_kjk[:,k]
   
           @. dy_p += c_jk[k]*beta_jk[k]*exp((beta_jk[k]*t)+phi_jk[k])*k_jk
   
       end
   
       # Returns the response
       return nothing
   
   end
   