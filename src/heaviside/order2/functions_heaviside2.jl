

"""
Driver to evaluate the permanent solution at t 

Inputs:

t           -> time

``sol_j``  -> cache computed by Process_heaviside

``load_data`` -> OrderedDict with loading information

CbF = Cb - F211 

M01 = M0^(-1)

M02 = M0^(-2)

M03 = M0^(-3)

M1 = F211^(-1)

M2 = F211^(-2)

M3 = F211^(-3)

M001 = Cb2F^(-1)

F211  -> constant matrix

m01m1 = M01*M1

m01m2 = M01*M2

m01m3 = M01*M3

m02m1 = M02*M1

m02m2 = M02*M2

m03m1 = M03*M1

outp         -> output vector (modified in place) 
"""
function y_permanent_heaviside2!(t::Float64,sol_j::AbstractMatrix,load_data::OrderedDict,CbF::AbstractMatrix,
                      M1::AbstractMatrix,
                      M2::AbstractMatrix, M3::AbstractMatrix,  M001,
                      F211::AbstractMatrix, m01m1::AbstractMatrix, m01m2::AbstractMatrix, m01m3::AbstractMatrix,
                      m02m1::AbstractMatrix, m02m2::AbstractMatrix, m03m1::AbstractMatrix,
                      outp::Vector{Ts})  where Ts

    
    # Set outp to zero
    fill!(outp,zero(Ts))                  

    # Number of cached entries
    ncol = size(sol_j,2)

    # 
    # Compute Eq. 179
    #
    for j in keys(load_data)

        # Data
        data = load_data[j]

        # Number of terms
        nk = Int(length(data)/4)

        # Counter (for j)
        count = 1

        # Loop over k
        for k=0:nk-1 

            # Polynomial Coefficients
            c_jk0 = data[4*k+1]  
            c_jk1 = data[4*k+2]  
            c_jk2 = data[4*k+3]  

            # offset
            t_jk = data[4*k+4]  
            
            # Heaviside
            H_jk = Heaviside(t,t_jk)

            # Just evaluate if H is not null
            if H_jk > 0.0

                # Common term
                t1 = t_jk - t

                # T1
                T1 = (c_jk2*t^2)*m01m1

                # T2
                T2 = t*(c_jk1*m01m1 .- 2*c_jk2*m01m2 .- 2*c_jk2*m02m1) 

                # T3 
                T3 = 2*c_jk2*m03m1

                # T4
                T4 = -c_jk1*m02m1 .+ 2*c_jk2*m02m2 
 
                # T5
                T5 = 2*c_jk2*m01m3 .- c_jk1*m01m2 .+ c_jk0*m01m1

                # T6
                T6A = -( c_jk2*t_jk^2 .+ c_jk1*t_jk .+ c_jk0  )*M1
                T6B =  (2*c_jk2*t_jk .+ c_jk1 )*M2
                T6C =  -2*c_jk2*M3
                T6 = exp(F211*t1)*(M001\( T6A .+ T6B .+ T6C))

                # T7
                T7A1 = m01m1*(-c_jk2*t_jk^2  -c_jk1*t_jk -c_jk0)   
                T7A2 = m01m2*(2*c_jk2*t_jk + c_jk1) 
                T7A3 = m01m3*(-2*c_jk2)
                T7A  = T7A1 .+ T7A2 .+ T7A3

                T7B1 =  m02m1*(2*c_jk2*t_jk + c_jk1)
                T7B2 =  m02m2*(-2*c_jk2)
                T7B  =  T7B1 .+ T7B2

                T7C = -T3 #-2*c_jk2*m03m1

                T7D1 = -M1*(c_jk2*t_jk^2 + c_jk1*t_jk + c_jk0)
                T7D2 =  M2*(2*c_jk2*t_jk + c_jk1)
                T7D3 = T6C #-2*c_jk2*M3
                T7D  = M001\-(T7D1 .+ T7D2 .+ T7D3)

                # Final T7
                T7 = exp(CbF*t1)*(T7A .+ T7B .+ T7C .+ T7D)

                # Add and use the cache
                outp .= outp .+ H_jk*( T1 .+ T2 .+ T3 .+ T4 .+ T5 .+ T6 .+ T7)*sol_j[:,count]

            end # if H_jk > 0.0
            
        end #k

        # Update count
        count += 1

    end #j

end

"""
Evaluate the permanent solution at t 

Inputs:

t           -> time

``sol_j``  -> cache computed by Process_heaviside

``load_data`` -> OrderedDict with loading information

CbF = Cb - F211 

M01 = M0^(-1)

M02 = M0^(-2)

M03 = M0^(-3)

M1 = F211^(-1)

M2 = F211^(-2)

M3 = F211^(-3)

M001 = Cb2F^(-1)

F211  -> constant matrix

m01m1 = M01*M1

m01m2 = M01*M2

m01m3 = M01*M3

m02m1 = M02*M1

m02m2 = M02*M2

m03m1 = M03*M1

Output: 

outp         -> output vector
"""
function y_permanent_heaviside2(t::Float64,sol_j::AbstractMatrix,load_data::OrderedDict,CbF::AbstractMatrix,
                    M1::AbstractMatrix,
                    M2::AbstractMatrix, M3::AbstractMatrix,  M001,
                    F211::AbstractMatrix, m01m1::AbstractMatrix, m01m2::AbstractMatrix, m01m3::AbstractMatrix,
                    m02m1::AbstractMatrix, m02m2::AbstractMatrix, m03m1::AbstractMatrix)


    # Number of DOFs
    ngls = size(sol_j,1)

    # Alocate the output vector
    outp = Vector{ComplexF64}(undef,ngls)

    # Call the driver
    y_permanent_heaviside2!(t,sol_j,load_data,CbF,M1,M2,M3,M001,F211,m01m1,m01m2,m01m3,
                           m02m1,m02m2,m03m1,outp)

    # return the output
    return outp

end
