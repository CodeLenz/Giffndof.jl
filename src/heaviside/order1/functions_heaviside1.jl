
"""
Driver to evaluate the permanent solution at t for Heavisides multiplied
    by first order polynomials

Inputs:

t           -> time

``sol_j``  -> cache computed by Process_heaviside

``load_data`` -> OrderedDict with loading information

CbF = Cb - F211 

M01 = M0^(-1)

M1 = F211^(-1)

M2 = F211^(-2)

M001 = Cb2F^(-1)

F211  -> constant matrix

m01m1 = M01*M1

m01m2 = M01*M2

m02m1 = M02*M1

outp         -> output vector (modified in place) 
"""
function y_permanent_heaviside1!(t::Float64,sol_j::AbstractMatrix,load_data::OrderedDict,CbF::AbstractMatrix,
                      M01::AbstractMatrix, M1::AbstractMatrix,
                      M2::AbstractMatrix, M001::AbstractMatrix,
                      F211::AbstractMatrix, m01m1::AbstractMatrix, m01m2::AbstractMatrix,
                      m02m1::AbstractMatrix, outp::Vector{Ts})  where Ts

    
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
        nk = Int(length(data)/3)

        # Counter (for j)
        count = 1

        # Loop over k
        for k=0:nk-1 

            # Polynomial Coefficients
            c_jk0 = data[3*k+1]  
            c_jk1 = data[3*k+2]  
         
            # offset
            t_jk = data[3*k+3]  
            
            # Heaviside
            H_jk = Heaviside(t,t_jk)

            # Just evaluate if H is not null
            if H_jk > 0.0

                # Common term
                t1 = t_jk - t

                # T2
                T2 = t*(c_jk1*m01m1) 

                # T4
                T4 = -c_jk1*m02m1
 
                # T5
                T5 = -c_jk1*m01m2 + c_jk0*m01m1

                # T6
                T6A = -(c_jk1*t_jk + c_jk0)*M1
                T6B =  (c_jk1)*M2
                T6 = exp(F211*t1)*M001*( T6A .+ T6B)

                # T7
                T7A1 = m01m1*(-c_jk1*t_jk -c_jk0)   
                T7A2 = m01m2*(c_jk1) 
                T7A  = T7A1 .+ T7A2

                T7B =  m02m1*(c_jk1)
             
                T7D1 = -M1*(c_jk1*t_jk + c_jk0)
                T7D2 =  M2*(c_jk1)
                T7D  = -M001*(T7D1 .+ T7D2)

                # Final T7
                T7 = exp(CbF*t1)*(T7A .+ T7B .+ T7D)

                # Add and use the cache
                outp .= outp .+ H_jk*( T2 .+ T4 .+ T5 .+ T6 .+ T7)*sol_j[:,count]

            end # if H_jk > 0.0
            
        end #k

        # Update count
        count += 1

    end #j

end

"""
Evaluate the permanent solution at t for Heavisides multiplied
    by first order polynomials

Inputs:

t           -> time

``sol_j``  -> cache computed by Process_heaviside

``load_data`` -> OrderedDict with loading information

CbF = Cb - F211 

M01 = M0^(-1)

M1 = F211^(-1)

M2 = F211^(-2)

M001 = Cb2F^(-1)

F211  -> constant matrix

m01m1 = M01*M1

m01m2 = M01*M2

m02m1 = M02*M1

Output: 

outp         -> output vector
"""
function y_permanent_heaviside1(t::Float64,sol_j::AbstractMatrix,load_data::OrderedDict,CbF::AbstractMatrix,
                    M01::AbstractMatrix, M1::AbstractMatrix,
                    M2::AbstractMatrix, M001::AbstractMatrix,
                    F211::AbstractMatrix, m01m1::AbstractMatrix, m01m2::AbstractMatrix, 
                    m02m1::AbstractMatrix)


    # Number of DOFs
    ngls = size(sol_j,1)

    # Alocate the output vector
    outp = Vector{ComplexF64}(undef,ngls)

    # Call the driver
    y_permanent_heaviside1!(t,sol_j,load_data,CbF,M01,M1,M2,M001,F211,m01m1,m01m2,
                           m02m1,outp)

    # return the output
    return outp

end

