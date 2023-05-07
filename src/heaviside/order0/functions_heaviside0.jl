
"""
Driver to evaluate the permanent solution at t for Heavisides multiplied
    by zero order polynomials

Inputs:

t           -> time

``sol_j``  -> cache computed by Process_heaviside

``load_data`` -> OrderedDict with loading information

CbF = Cb - F211 

M01 = M0^(-1)

M1 = F211^(-1)

M001 = Cb2F^(-1)

F211  -> constant matrix

m01m1 = M01*M1

outp         -> output vector (modified in place) 
"""
function y_permanent_heaviside0!(t::Float64,sol_j::AbstractMatrix,load_data::OrderedDict,CbF::AbstractMatrix,
                      M1::AbstractMatrix,  M001,
                      F211::AbstractMatrix, m01m1::AbstractMatrix, 
                      outp::Vector{Ts})  where Ts

    
    # Set outp to zero
    fill!(outp,zero(Ts))                  

    # Number of cached entries
    ncol = size(sol_j,2)

    # 
    # Compute Eq. 
    #
    for j in keys(load_data)

        # Data
        data = load_data[j]

        # Number of terms
        nk = Int(length(data)/2)

        # Counter (for j)
        count = 1

        # Loop over k
        for k=0:nk-1 

            # Polynomial Coefficients
            c_jk0 = data[2*k+1]  
         
            # offset
            t_jk = data[2*k+2]  
            
            # Heaviside
            H_jk = Heaviside(t,t_jk)

            # Just evaluate if H is not null
            if H_jk > 0.0

                # Common term
                t1 = t_jk - t
 
                # T5
                T5 = c_jk0*m01m1

                # T6
                T6A = -c_jk0*M1
                CT = M001\T6A
                T6  = exp(F211*t1)*CT

                # T7
                #T7A = -T5 #m01m1*(-c_jk0)   <------------------
                
                #T7D1 = T6A #-M1*(c_jk0)  <----------
                T7D  = -CT

                # Final T7
                T7 = exp(CbF*t1)*(T7D.-T5)

                # Add and use the cache
                outp .= outp .+ H_jk*(T5 .+ T6 .+ T7)*sol_j[:,count]

            end # if H_jk > 0.0
            
        end #k

        # Update count
        count += 1

    end #j

end

"""
Evaluate the permanent solution at t for Heavisides multiplied
    by zero order polynomials

Inputs:

t           -> time

``sol_j``  -> cache computed by Process_heaviside

``load_data`` -> OrderedDict with loading information

CbF = Cb - F211 

M01 = M0^(-1)

M1 = F211^(-1)

M001 = Cb2F^(-1)

F211  -> constant matrix

m01m1 = M01*M1

Output: 

outp         -> output vector
"""
function y_permanent_heaviside0(t::Float64,sol_j::AbstractMatrix,load_data::OrderedDict,CbF::AbstractMatrix,
                    M1::AbstractMatrix,
                    M001, F211::AbstractMatrix, m01m1::AbstractMatrix)


    # Number of DOFs
    ngls = size(sol_j,1)

    # Alocate the output vector
    outp = Vector{ComplexF64}(undef,ngls)

    # Call the driver
    y_permanent_heaviside0!(t,sol_j,load_data,CbF,M1,M001,F211,m01m1,outp)

    # return the output
    return outp

end
