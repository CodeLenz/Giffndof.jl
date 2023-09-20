function polynomial_particularResponse(times, n_times, data, response, tol_complex=1E-6)

        # Initializes a counter of accessed vectors v_klp

        counter_vklp = 1

        # Iterates through the k polynomial terms
        for i=1:data.nk
    
            # Access the degree
            k = data.degrees[i]
    
            # Access the time shift
            t_j = data.t_j[i]
    
            # Access the monomial coefficients
            c_jk = @view data.c_j[i]

            # Iterates through the summation in l
            for l=1:(k+1)
    
                # Iterates through the summation in p
                for p=1:(k-l+2)
    
                    # Access the vector v_k_l_p
                    v_klp = @view data.matrix_vklp[:,counter_vklp]

                    # Updates the counter
                    counter_vklp += 1
    
                    # Iterates through the time points
                    for m=1:n_times
    
                        @. response[:,m] += c_jk*((times[m]-t_j)^(k-l-p+2))*v_klp
    
                    end
    
                end 
    
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
    # tion due to polynomial excitation using the data in a dictionary
    function polynomial_derivativeParticular(t, data, dy_p)
    
        # Initializes a counter of accessed vectors v_klp

        counter_vklp = 1

        # Iterates through the k polynomial terms
        for i=1:data.nk
    
            # Access the degree
            k = data.degrees[i]
    
            # Access the time shift
            t_j = data.t_j[i]
    
            # Access the monomial coefficients
            c_jk = @view data.c_j[i]
    
            # Iterates through the summation in l
            for l=1:(k+1)
    
                # Iterates through the summation in p
                for p=1:(k-l+2)
    
                    # Tests whether this monomial in the solution is non-
                    # constant and, if so, derivates it in conventional fa-
                    # shion
                    if (k-l-p+2)>0
    
                        # Access the vector v_k_l_p
                        v_klp = @view data.matrix_vklp[:,counter_vklp]
    
                        # Iterates through the time points
                        @. dy_p += c_jk*(k-l-p+2)*((t-t_j)^(k-l-p+1))*v_klp
    
                    end

                    # Updates the counter out of the if statement to not
                    # lose pace even though the vector v_klp is not used

                    counter_vklp += 1
    
                end 
    
            end
    
        end 
    
        return nothing

    end
    