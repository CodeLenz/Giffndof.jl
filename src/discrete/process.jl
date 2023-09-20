#=

Gera dois dicionários contendo as informações pós-processadas sobre os carregamentos

data -> coeficientes e soluções dos sistemas de equações para cada tipo de solução

excitations -> rotinas para as soluções particulares

=#

function Process(dimen, M, C, K, load_description,F11, C_bar, K_bar,
 input_vector, norm_conjugacy, expF11_delta, expCF_delta, CbF, tol)
   
       ##
       ## Exponential
       ##
   
       # Initializes the vectors for the coefficients of the exponential terms
       exp_c_jk = ComplexF64[]
       exp_beta_jk = ComplexF64[]
       exp_phi_jk = ComplexF64[]

       # Initializes a matrix for the vectors k_jk of the exponential
       matrix_kjk = zeros(ComplexF64, dimen)

       ##
       ## Polynomial
       ##

       # Initializes the vectors of coefficients, time shifts and powers
       pol_c_j = Float64[]
       pol_t_j = Float64[]
       pol_degree = Int64[]

       # Initializes a matrix for the vectors v_klp of the polynomial

       matrix_vklp = zeros(ComplexF64, dimen)

       ##
       ## Dirac
       ##

       # Initializes a matrix for coefficients and times of the Dirac's 
       # deltas to guarantee that they will be placed in a crescent order
       # of time application
       diracs_timesNCoeffs = [0.0; 0.0]

       dir_vj = zeros(ComplexF64,dimen)

       ##
       ## Miscelaneous
       ##

       # Initializes counters of the number of terms for each of the possible excitations, excitation_nk
       exp_nk = 0
       pol_nk = 0
       dir_nk = 0
   
       # Initializes a dictionary of flags and functions to signal which
       # excitation types has been selected. The first column has the func-
       # tion for the particular response; the second column, the function
       # for the derivative of the particular response; the third column,
       # the flag itself. The elements are ordered for exponential, polyno-
       # mial and Dirac's delta
      
       # Exponential
       flag_exp = false 
   
       # Polynomial
       flag_pol = false 
   
       # Dirac's delta
       flag_dir = false 

       # Iterates through the loads
       for i=1:size(load_description,1)
   
           # If it is an exponential
           if load_description[i,1]=="exponential"
   
               # Updates the flag for excitation
               # excitations["exp"][3] = true
               flag_exp = true
   
               # Allocates the coefficients
               push!(exp_c_jk,load_description[i,2])
               push!(exp_beta_jk,load_description[i,3])
               push!(exp_phi_jk,load_description[i,4])
   
               # Updates the exponential_nk counter
               exp_nk += 1
   
               # angular frequency
               angular = load_description[i,3]

               # Calculates k_jk and alocates it into the matrix
               matrix_kjk = [matrix_kjk ((K+(angular*C)+(angular^2)*M)\input_vector)]
   
           # Otherwise, if it is a sine
           elseif load_description[i,1]=="sine" || load_description[i,1]=="cosine"
               
               # Updates the flag for excitation
               # excitations["exp"][3] = true
               flag_exp = true
   
               # Allocates the coefficients
               if load_description[i,1]=="sine"
                    push!(exp_c_jk, load_description[i,2]*0.5im)
                    push!(exp_beta_jk, -2*pi*load_description[i,3]*1im)
                    push!(exp_phi_jk, -(pi/180)*load_description[i,4]*1im)
                    push!(exp_c_jk, (-load_description[i,2]*0.5im))
                    push!(exp_beta_jk, 2*pi*load_description[i,3]*1im)
                    push!(exp_phi_jk, (pi/180)*load_description[i,4]*1im)
               else
                    push!(exp_c_jk, load_description[i,2]*0.5)
                    push!(exp_beta_jk, 2*pi*load_description[i,3]*1im)
                    push!(exp_phi_jk, (pi/180)*load_description[i,4]*1im)
                    push!(exp_c_jk, load_description[i,2]*0.5)
                    push!(exp_beta_jk, -2*pi*load_description[i,3]*1im)
                    push!(exp_phi_jk, -(pi/180)*load_description[i,4]*1im)     
               end
                    
               # Updates the exponential_nk counter
               exp_nk += 1
   
               # Angular frequency
               angular = exp_beta_jk[end-1]

               # Calculates k_jk and alocates it into the matrix
               matrix_kjk = [matrix_kjk ((K.+(angular*C).+(angular^2)*M)\input_vector)]
   
               # Updates the exponential_nk counter
               exp_nk += 1
   
               # Calculates k_jk and alocates it into the matrix
               matrix_kjk = [matrix_kjk conj(matrix_kjk[:,end])]
   
            elseif load_description[i,1]=="polynomial"

                # Updates the flag for excitation    
                # excitations["pol"][3] = true
                flag_pol = true
    
                # Updates the polynomial_nk counter
                pol_nk += 1
    
                # Creates the vectors v_k_l_p, that are equal to [M*(F11^(l-p))*(K_bar^p)]\e_j
    
                # Gets the degree of the monomial term and saves it into the
                # dictionary
                k = load_description[i,4]
   
                # Tests if k is an integer
                if (k-floor(k))>0
    
                    throw("The exponent "*string(k)*" is not an integer, which is not allowed\n")
    
                else
    
                    # Turns it into an integer
                    k = Int(k)
    
                end

                # Adds the power to the vector of degrees
                push!(pol_degree, k)
    
                # Saves the coefficient of the polynomial and the shift
                push!(pol_c_j, load_description[i,2])
                push!(pol_t_j, load_description[i,3])
    
                # Iterates through the summation in l
                for l=1:(k+1)
    
                    # Calculates the term with factorials
                    factorial_l = ((-1)^(l+1))*(factorial(k)/factorial(k-l+1))
    
                    # Iterates through the summation in p
                    for p=1:(k-l+2)
    
                        # Calculates the term with factorials
                        factorial_p = ((-1)^(p+1))*(factorial(k-l+1)/factorial(k-l-p+2))
    
                        # Calculates the vector v_k_l_p
                        v_klp = (M*(F11^(l-p))*(K_bar^p))\input_vector
    
                        # Multiplates v_k_l_p by the factorial terms and al-
                        # locates it into the matrix
                        matrix_vklp = [matrix_vklp (v_klp*factorial_l*factorial_p)]
    
                    end #p
    
                end #l
    
            elseif load_description[i,1]=="dirac"

                # Updates the flag for excitation
                flag_dir = true
                #excitations["dir"][3] = true
    
                # Updates the dirac_nk counter
                dir_nk += 1
    
                # Saves the coefficient and the time of application
                diracs_timesNCoeffs = [diracs_timesNCoeffs (load_description[i,2:3])]
    
           # Otherwise, there is no option but throw error
           else
               throw(("Option of load "*load_description[i,1]*" is not available\n"))
           end
   
       end

    
    # Saves the struct for the exponential function
    if flag_exp 

        data_exponential = exponential_dataStruct(exp_nk, exp_c_jk,
         exp_beta_jk, exp_phi_jk, matrix_kjk[:,2:end])

    else 

        data_exponential = exponential_dataStruct(exp_nk, [0.0],
        [0.0], [0.0], [0.0 0.0])

    end

    # Saves the struct for the polynomial function
    if flag_pol 

        data_polynomial = polynomial_dataStruct(pol_nk, pol_c_j, pol_t_j,pol_degree, matrix_vklp[:,2:end])

    else

        data_polynomial = polynomial_dataStruct(pol_nk, [0.0;0.0], [0.0;0.0], [0;0], [0.0 0.0])

    end
   
    # Saves the input vector for Dirac's delta and the exponentials
    if flag_dir 

        # Calculates the vector v_j
        dir_vj .= (C.-(2*M*F11))\input_vector

        # Takes out the first column of the diracs_timesNCoeffs
        diracs_timesNCoeffs = diracs_timesNCoeffs[:,2:end]

        # Orders the matrix of coefficients and times in crescent order
        # of time applications
        diracs_timesNCoeffs = diracs_timesNCoeffs[:,sortperm(diracs_timesNCoeffs[2,:])]

        # Adds the flag for conjugacy of Cb-F11 and F11
        if norm_conjugacy<tol

            # Creates the struct
            data_dirac = dirac_dataStructConjugate(dimen, dir_nk,
             diracs_timesNCoeffs[1,:], diracs_timesNCoeffs[2,:], dir_vj,
             true, F11, expF11_delta)

        else

            # Creates the struct
            data_dirac = dirac_dataStructNonConjugate(dimen, dir_nk,
             diracs_timesNCoeffs[1,:], diracs_timesNCoeffs[2,:], dir_vj,
             false, F11, CbF, expF11_delta, expCF_delta)

        end

    else

        data_dirac = dirac_dataStructConjugate(dimen, dir_nk,[0.0;0.0], [0.0;0.0], [0.0;0.0], true, [0.0 0.0], [0.0 0.0])


    end

    # Returns the dictionary  
    return data_exponential, data_polynomial, data_dirac, flag_exp, flag_pol, flag_dir
   
   end
   