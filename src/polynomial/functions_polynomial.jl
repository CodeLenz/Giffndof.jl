#
# Solution for a polynomial excitation, Eq 150.
#
#  ``g_j(t) = c_j0 (t-t_j)^0 + c_j1 (t-t_j)^1 + ... + c_j2k (t-t_j)^k``
#
# where j is the DOF and k=1..nk is the number of base terms
#


"""
Pre-process data to avoid some computations when evaluating the 
permanent solution yp(t). Those computations are used many times,
so we evaluate them and store into a cache.

Inputs:

 M -> Constant matrix

 Kb -> Constant matrix ``Kb = M\\K``

 F211 -> Constant matrix

 ``load_data`` -> Dictionary with key j (gl) and data ``[t_j c_j0 ...  c_jnk]``
 
Outputs:

 ``sol_jlp`` -> Vector of vectors contaning ``M Kb^p F^{l-p}_{2,1,2} \\ e_j``

"""
function Process_polynomial(M::AbstractMatrix{T}, Kb::AbstractMatrix{T}, F211::AbstractMatrix,
                      load_data::OrderedDict{Int64,Vector{Float64}}) where T


    # Number of DOFs
    ngls = size(M,1)

    # Cache
    sol_jklp = Vector{ComplexF64}[]

    # Pre allocate e_j
    e_j = zeros(ngls)

    #
    # Compute M Kb^p F^(l-p) \ e_j
    #
    cont = 1
    for j in keys(load_data)

        # Put 1.0 at j
        e_j[j] = 1.0

        # Coefficients c_jk 
        data = load_data[j]
        coefs = data[2:end]

        # Find nk 
        nk = length(coefs)

        # Loop in k
        for k=0:nk-1

            # Coefficient c_jk
            c_jk = coefs[k+1]

            # Avoid computations for null coefficients
            if c_jk != 0.0

                # Loop in l
                for l=1:k+1

                    # Loop in p 
                    for p=1:k-l+2

                            # Evaluate M Kb^p F211^(l-p) 
                            Kjklp = M*((F211^(l-p))*(Kb^p))

                            # Store in sol_jkl
                            push!(sol_jklp,Kjklp\e_j)

                    end #p   

                end #l

            end #else c_jk != 0.0

        end #k

        # Unset o e_j 
        e_j[j] = 0.0

    end #j

    # Return cache sol_jklp
    return sol_jklp

end

"""
Driver to evaluate the permanent solution at t 

Inputs:

t           -> time

``sol_jklp``  -> cache computed by Process_poly

``load_data`` -> OrderedDict with loading information

outp         -> output vector (modified in place) 
"""
function y_permanent_polynomial!(t::Float64,sol_jklp::Vector{Vector{T}},
                           load_data::OrderedDict{Int64,Vector{Float64}},
                           outp::Vector{T}) where T

    # Make sure outp is zero
    fill!(outp,zero(T))

    #
    # Compute Eq. 161
    #
    cont = 1
    for j in keys(load_data)

        # Get data
        data = load_data[j]

        # get time t_j
        tj = data[1]

        # Coefficients
        coefs = data[2:end]

        # Number of terms nk
        nk = length(coefs)

        # Loop in k
        for k=0:nk-1

            # Coefficient c_jk
            c_jk = coefs[k+1]

            # Avoid computation for null c_jk
            if c_jk != 0.0

                # Loop in l
                for l=1:k+1

                    # Coefficient for l
                    cl_1 = (-1)^(l+1)
                    cl_2 = factorial(k)/factorial(k-l+1)
                    cl = cl_1 * cl_2

                    # Loop in p 
                    for p=1:k-l+2

                        # Coefficient
                        cp_1 = (-1)^(p+1)
                        cp_2 = factorial(k-l+1)/factorial(k-l-p+2)
                        cp_3 = (t-tj)^(k-l-p+2)
                        cp = cp_1 * cp_2 * cp_3

                        # Use the cache and add to the solution
                        outp .= outp .+ c_jk*cl*cp*sol_jklp[cont]

                        # update cont
                        cont += 1
                        
                    end #p   

                end #l

            end # if c_jk != 0.0

        end #k

    end #j
end

"""
Evaluate the permanent solution at t 

Inputs:

t           -> time

``sol_jklp``  -> cache computed by Process_poly

``load_data`` -> OrderedDict with loading information

Outputs:

outp         -> output vector
"""
function y_permanent_polynomial(t::Float64,sol_jklp::Vector{Vector{T}},
                          load_data::OrderedDict{Int64,Vector{Float64}}) where T

    # Number of DOFs
    ngls = length(sol_jklp[1])

    # Alocate output vector
    outp = Vector{T}(undef,ngls)

    # Call the driver
    y_permanent_polynomial!(t,sol_jklp,load_data,outp)

    # return 
    return outp

end


"""
Driver to evaluate the derivative of the permanent solution at t 

Inputs:

t           -> time

``sol_jklp``  -> cache computed by Process_poly

``load_data`` -> OrderedDict with loading information

outp         -> output vector (modified in place) 
"""
function dy_permanent_polynomial!(t::Float64,sol_jklp::Vector{Vector{T}},
                            load_data::OrderedDict{Int64,Vector{Float64}},
                            outp::Vector{T}) where T

    # Make sure outp is zero
    fill!(outp,zero(T))

    # 
    # Compute Eq. 162
    #
    cont = 1
    for j in keys(load_data)

     # Get data
     data = load_data[j]

     # get time t_j
     tj = data[1]

     # Coefficients
     coefs = data[2:end]

     # Number of terms nk
     nk = length(coefs)

     # Loop in k
     for k=0:nk-1

         # Coefficient c_jk
         c_jk = coefs[k+1]

         # Avoid computations for null c_jk
         if c_jk != 0.0

            # Loop in l
            for l=1:k+1

                # Coefficient for l
                cl_1 = (-1)^(l+1)
                cl_2 = factorial(k)/factorial(k-l+1)
                cl = cl_1 * cl_2

                # Loop in p 
                for p=1:k-l+2

                    # The main difference for yp(t)
                    cp_4 = k-l-p+2 

                    # There is no point in evaluating it if cp_4 == 0
                    if cp_4 != 0.0

                        # Coefficient
                        cp_1 = (-1)^(p+1)
                        cp_2 = factorial(k-l+1)/factorial(k-l-p+2)
                        cp_3 = (t-tj)^(k-l-p+1)

                        # Final coefficient
                        cp = cp_1 * cp_2 * cp_3 * cp_4

                        # Now we can add and use the cache
                        outp .= outp .+ c_jk*cl*cp*sol_jklp[cont]
                   
                    end # k-l-p+2

                    # update cont
                    cont += 1

                end #p   

            end #l

        end # if c_jk != 0.0

     end #k

    end #j

end

"""
Evaluate the derivative of the permanent solution at t 

Inputs:

t           -> time

``sol_jklp``  -> cache computed by Process_poly

``load_data`` -> OrderedDict with loading information

Outputs:

outp         -> output vector
"""
function dy_permanent_polynomial(t::Float64,sol_jklp::Vector{Vector{T}},
                           load_data::OrderedDict{Int64,Vector{Float64}}) where T

    # Number of DOFs
    ngls = length(sol_jklp[1])

    # Allocate the output
    outp = Vector{T}(undef,ngls)

    # Call the driver
    dy_permanent_polynomial!(t,sol_jklp,load_data,outp)

    # Return the solution at t
    return outp

end

