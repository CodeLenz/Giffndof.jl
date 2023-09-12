#=
 Solution procedure to systems of second order ODEs 

 M a(t) + C v(t) + K y(t) = F(t)

 Where M, C and K are n x n matrices 

 y(t) is n x 1 -> dependent variable

 v(t) is n x 1 -> the first derivative of y(t)

 a(t) is n x 1 -> the second derivative of y(t)

 and f(t) is n x 1 with excitations in the form

  f_j(t) = sum_k  c_jk exp( beta_jk t + im*phi_jk)

 where j is the DOF, c_jk is a complex amplitude, 
 beta_jk = omega_jk * i is a complex angular frequency and phi_jk the phase.

 Loading information is given by using a dictionary:

 load_data -> dictionary with key j (gl) and 
 data [c_j1; w_j1; phi_j1; c_j2; w_j2; phi_j2; ... c_jnk; w_jnk; phi_jnk]
=#

"""
  Pre-process data to avoid some computations when evaluating the 
  permanent solution yp(t). Those computations are used many times,
  so we evaluate them and store into three caches: 

   ``c_jk*(KD_jk \\ e_j) -> sol_jk`` 

 
 ``im*w_jk -> beta_jk`` 

 and

 ``im*phi_jk -> cphi_jk`` 


 Inputs:

 M         -> Constant matrix

 C         -> Constant matrix

 K         -> Constant matrix

 ``load_data`` -> dictionary with key j (gl) and data ``[c_j1; w_j1; phi_j1; ... c_jnk; w_jnk; phi_jnk]``

 Outputs

 ``sol_jk``  -> complex matrix

 ``beta_jk`` -> complex vector

 ``cphi_jk`` -> complex vector


"""
function Process_exponential(M::AbstractMatrix{T}, C::AbstractMatrix{T},
                     K::AbstractMatrix{T}, 
                     load_data::Dict{Int64,Vector{ComplexF64}}) where T


    # Basic assertions
    @assert all(size(M).==size(C).==size(K)) "Process_exponential:: Coeficient matrices do not have the same size"
    @assert !isempty(load_data) "Process_exponential:: There is no loading data"

    # Number of DOFs
    ngls = size(M,1)

    # Loop over the dictionary to find the total number of data to be cached
    ncol = 0
    for j in keys(load_data)
        ncol += Int(length(load_data[j])/3)
    end

    # Allocate caches  
    sol_jk  = zeros(ComplexF64,ngls,ncol)
    beta_jk = zeros(ComplexF64,ncol)
    cphi_jk = zeros(ComplexF64,ncol)

    # Also allocate vetor e_j 
    e_j = zeros(ngls)

    # For each DOF j in load_data, evaluate c_jk*(KD_jk\e_j)
    # and im*w_jk 
    cont = 1
    for j in keys(load_data)

        # Data
        data = load_data[j]

        # Number of data for j is 
        nk = Int(length(data)/3)

        # Put 1.0 in DOF j
        e_j[j] = 1.0

        # Loop over k 
        for k=1:nk

            # Amplitude
            c_jk = data[3*(k-1)+1]

            # Angular frequency
            w_jk = real(data[3*(k-1)+2])

            # Phase
            phi_jk = real(data[3*(k-1)+3])

            # Build KD_jk
            KD_jk = K + im*w_jk*C - M*(w_jk)^2 

            # Solve for  K_jk\e_j and multiply by c_jk
            sol_jk[:,cont] .= c_jk*(KD_jk\e_j)
           
            # Store im*w_jk
            beta_jk[cont] = im*w_jk 

            # Store im*phi_jk
            cphi_jk[cont] = im*phi_jk 

            # Increment the counter
            cont += 1

        end #k
        
        # Unset e_j 
        e_j[j] = 0.0

    end #j

    # Return the caches
    return sol_jk, beta_jk, cphi_jk

end


"""
Driver to evaluate the permanent solution at t 

Inputs 

t           -> time

``sol_jk``  -> cache computed by Process_exp

``beta_jk`` -> cache computed by Process_exp

``cphi_jk`` -> cache computed by Process_exp

outp         -> output vector (modified in place) 
"""
function y_permanent_exponential!(t::Float64,sol_jk::AbstractMatrix{T},beta_jk::Vector{T},cphi_jk::Vector{T},
                                  outp::Vector{T}) where T

    # Number of informations stored in the caches
    ncol = size(sol_jk,2)

    # Set the output to zero
    fill!(outp,zero(T))

    # Loop - Equation  
    for jk=1:ncol
        outp .= outp .+ exp(beta_jk[jk]*t + cphi_jk[jk])*sol_jk[:,jk]
    end

end


"""
Evaluate the permanent solution at t 

Inputs 

t           -> time

``sol_jk``  -> cache computed by Process_exp

``beta_jk`` -> cache computed by Process_exp

``cphi_jk`` -> cache computed by Process_exp

Output

otp         -> output vector
"""
function y_permanent_exponential(t::Float64,sol_jk::AbstractMatrix{T},beta_jk::Vector{T},
                                 cphi_jk::Vector{T} ) where T

    # Number of DOFs
    ngls = size(sol_jk,1)

    # Create the output vector
    outp = Vector{T}(undef,ngls)

    # Call the driver
    y_permanent_exponential!(t,sol_jk,beta_jk,cphi_jk,outp)
    
    # Return the output
    return outp

end


"""
Driver to evaluate the derivative of the permanent solution at t 

Inputs 

t           -> time

``sol_jk``  -> cache computed by Process_exp

``beta_jk`` -> cache computed by Process_exp

``cphi_jk_jk`` -> cache computed by Process_exp

otp         -> output vector (modified in place) 
"""
function dy_permanent_exponential!(t::Float64,sol_jk::AbstractMatrix{T},beta_jk::Vector{T}, 
                                  cphi_jk::Vector{T},outp::Vector{T}) where T

    # Number of informations stored in the caches
    ncol = size(sol_jk,2)

    # Make sure outp is zero
    fill!(outp,zero(T))

    # Loop - Equation
    for jk=1:ncol
        outp .= outp .+ beta_jk[jk]*exp(beta_jk[jk]*t + cphi_jk[jk])*sol_jk[:,jk]
    end

end



"""
Evaluate the derivative of the permanent solution at t 

Inputs 

t           -> time

``sol_jk``  -> cache computed by Process_exp

``beta_jk`` -> cache computed by Process_exp

Output

otp         -> output vector
"""
function dy_permanent_exponential(t::Float64,sol_jk::AbstractMatrix{T},beta_jk::Vector{T},
                                  cphi_jk::Vector{T}) where T

    # Number of informations stored in the caches
    ngls = size(sol_jk,1)

    # Create the output vector
    outp = Vector{T}(undef,ngls)

    # Call the driver
    dy_permanent_exponential!(t,sol_jk,beta_jk,cphi_jk,outp)

    # Return the output
    return outp

end

