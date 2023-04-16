"""
 Evaluate the coefficients a used to compute coefficients c.

 Continuous version (load is a function of time)

 Inputs:

 g is the reference (known) function 

 tl is the initial time at the interval

 tL is the final time at the interval

 Outputs:

 a0 coefficient 
 a1 coefficient 

"""
function Evaluate_coefs_a(g::Function, tl::Float64, tL::Float64)

    # Testing purposes
    tL > tl || error("Evaluate_coefs_a:: tl must be smaller than tL")

    # Values of g(t) at the extremes of the interval
    gl = g(tl)
    gL = g(tL)

    # Evaluate dt for this interval
    dt = tL-tl

    # Coeficient a1
    a1 = (gL-gl)/dt

    # Integral of g in the interval. 
    # Evalute using quadgk package.
    integral = quadgk(g,tl,tL)[1] 

    # Coeficient a0
    a0 = integral/dt - a1*(tL^2 - tl^2)/(2*dt)

    # Retorn both coefficients
    return a0, a1

end


"""
 Evaluate the coefficients a used to compute coefficients c.

 Discrete version (load is given as a vector)

 Inputs:

 interval is the initial position to the interval
 
 g is a vector with discrete values

 Ts is a vector of discrite times or a StepRange

 Outputs:

 a0 coefficient 
 a1 coefficient 

"""
function Evaluate_coefs_a(interval::Int64,g::Vector{T}, Ts::T1) where {T,T1}

    # Recover the times
    tl = Ts[interval]
    tL = Ts[interval+1]

    # Values of g at the extremes of the interval
    gl = g[interval]
    gL = g[interval+1]

    # Evaluate dt for this interval
    dt = tL-tl

    # Coeficient a1
    a1 = (gL-gl)/dt

    # Integral of g in the interval. 
    integral = dt*(gl + 0.5*(gL-gl)) 

    # Coeficient a0
    a0 = integral/dt - a1*(tL^2 - tl^2)/(2*dt)

    # Retorn both coefficients
    return a0, a1

end


"""
 Evaluate the coefficients c. Continuous version.

 Inputs:

 g is the reference (known) function 

 Ts a list or StepRange with discrete times
 
 Outputs:

 cj0 coefficient 
 cj1 coefficient 
 
"""
function Evaluate_coefs_c(g::Function, Ts::T) where T

   # Number of intervals 
   n = length(Ts)-1

   # Allocate both output vectors
   cj0 = zeros(n)
   cj1 = zeros(n)

   # Intial time
   tl = Ts[1]

   # Loop over the intervals
   for interval=1:n

       # Final time at this interval
       tL = Ts[interval+1]
       
       # Evaluate coefs a for this interval
       a0, a1 = Evaluate_coefs_a(g,tl,tL)

       # Evaluate coefs c
       cj0[interval] = a0 - sum(cj0[1:interval-1])
       cj1[interval] = a1 - sum(cj1[1:interval-1])

       # Update tl
       tl = tL

   end # interval

   # Return the coefficients
   cj0, cj1

end

"""
 Evaluate the coefficients c. Continuous version.

 Inputs:

 g is a vector with discrete values

 Ts is a vector of discrite times or a StepRange

 Outputs:

 cj0 coefficient 
 cj1 coefficient 
 
"""
function Evaluate_coefs_c(g::Vector{T}, Ts::T1) where {T,T1}

    # Basic asseertion
    length(g)==length(Ts) || error("Evaluate_coefs_c:: vectors g and Ts must have the same length")

   # Number of intervals 
   n = length(Ts)-1

   # Allocate both output vectors
   cj0 = zeros(n)
   cj1 = zeros(n)

   # Loop over the intervals
   for interval=1:n
 
       # Evaluate coefs a for this interval
       a0, a1 = Evaluate_coefs_a(interval,g,Ts)

       # Evaluate coefs c
       cj0[interval] = a0 - sum(cj0[1:interval-1])
       cj1[interval] = a1 - sum(cj1[1:interval-1])

   end # interval

   # Return the coefficients
   cj0, cj1

end


#
# Evaluate all coefficients cj0 and cj1 
# and store in a OrderedDict
#
function Generate_Dict_c(load_data::OrderedDict{Int64,T},Ts::T1) where {T,T1}

    # Create the dictionary with c's for every DOF j informed in load_data
    dict_c = OrderedDict{Int64,Matrix{Float64}}()

    # Loop over load_data
    for j in keys(load_data)

        # Evaluate the Coefficients c
        cj0, cj1 = Evaluate_coefs_c(load_data[j], Ts)

        # Store in the new dictionary
        dict_c[j] = [cj0 cj1]
        
    end #j

    # Return the dictionary with the c's for every DOF
    return dict_c

end


#
# Build the approximation for g(t) at t
#
"""
 Evaluate the approximation for a given function g(t) at t. 
 This approximation is built by using the Heaviside Series with linear polynomials
 at each interval. The purpose of this subroutine is for verification purposes only,
 since the evaluation of the response y⁽¹⁾(t) is performed by using the coefficients
 c only.

 Inputs:

 t is the time to evaluate the approximation

 cj0 is a vector of coefficients (obtained by using Evaluate_coefs_c)

 cj1 is a vector of coefficients (obtained by using Evaluate_coefs_c)

 Ts is a vector or StepRange with the discrete times used to evaluate cj0 and cj1

 Output:

 the approximation at
 
"""
function Evaluate_gtilde(t::Float64, cj0::Vector{Float64}, cj1::Vector{Float64}, Ts::T) where T

   # Number of intervals in Ts
   n = length(Ts)-1

   # Basic test
   n>1 || error("Evaluate_gtilde:: you must inform at least two values in Ts")

   # Inicialize the output vector
   saida = 0.0

   # We do not need to evaluate all intervals. Build a logical mask
   # where t is smaller than Ts (0)
   mask = t .< Ts

   # Thus, the number of intervals to use in the loop is
   np = min(n+1 - sum(mask),n)

   # Loop over the valid intervals
   for interval=1:np

        # Initial time at the interval
        tl = Ts[interval]
          
        # Add the contribution for this interval
        saida += Heaviside(t,tl)*(cj0[interval] + cj1[interval]*t)
   
    end # interval
   
   # Return the approximation to g(t)
   return saida

end


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

m02m1 = M02*M1

m02m2 = M02*M2

outp         -> output vector (modified in place) 
"""
function y_permanent_HS1!(t::Float64,sol_j::AbstractMatrix,dict_c::OrderedDict{Int64,Matrix{Float64}},Ts::T,
                      CbF::AbstractMatrix,
                      M01::AbstractMatrix, M02::AbstractMatrix,  M1::AbstractMatrix,
                      M2::AbstractMatrix,   M001::AbstractMatrix,
                      F211::AbstractMatrix, m01m1::AbstractMatrix, m01m2::AbstractMatrix,
                      m02m1::AbstractMatrix, m02m2::AbstractMatrix,
                      outp::Vector{Tp})  where {T,Tp}

    
    # Set outp to zero
    fill!(outp,zero(Tp))                  

    # Number of cached entries
    ncol = size(sol_j,2)

    # We do not need to evaluate all intervals. Build a logical mask
    # where t is smaller than Ts (0)
    mask = t .< Ts

    # Number of intervals
    n = length(Ts) - 1

    # Thus, the number of intervals to use in the loop is
    np = min(n+1 - sum(mask),n)

    # 
    # Compute Eq. 219
    #
    for j in keys(dict_c)

        # Counter (for j)
        count = 1

        # Recover coefficients c for this j
        cj0 = @view dict_c[j][:,1]
        cj1 = @view dict_c[j][:,2]
   
        # Loop over k just for valid intervals in Ts
        for k=1:np

            # Polynomial Coefficients
            c_jk0 = cj0[k]
            c_jk1 = cj1[k]
            
            # offset (initial time of this interval)
            t_jk = Ts[k]
            
            # Common term
            t1 = t_jk - t

            # T2
            T2 = t*(c_jk1*m01m1) 

            # T4
            T4 = -c_jk1*m02m1
 
            # T5
            T5 = -c_jk1*m01m2 .+ c_jk0*m01m1

            # T6
            T6A = -(c_jk1*t_jk + c_jk0)*M1
            T6B =  (c_jk1)*M2
            T6 = exp(F211*t1)*M001*(T6A .+ T6B)

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
            outp .= outp .+ (T2 .+ T4 .+ T5 .+ T6 .+ T7)*sol_j[:,count]

        end #k

        # Update count
        count += 1

    end #j

end


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

m02m1 = M02*M1

m02m2 = M02*M2

outp         -> output vector (modified in place) 
"""
function y_permanent_HS1(t::Float64,sol_j::AbstractMatrix,dict_c::OrderedDict{Int64,Matrix{Float64}},Ts::T,
                      CbF::AbstractMatrix,
                      M01::AbstractMatrix, M02::AbstractMatrix,  M1::AbstractMatrix,
                      M2::AbstractMatrix,   M001::AbstractMatrix,
                      F211::AbstractMatrix, m01m1::AbstractMatrix, m01m2::AbstractMatrix,
                      m02m1::AbstractMatrix, m02m2::AbstractMatrix)  where {T}



    # Number of DOFs
    ngls = size(sol_j,1)

    # Alocate the output vector
    outp = Vector{ComplexF64}(undef,ngls)

    # Call the driver
    y_permanent_HS1!(t,sol_j,dict_c,Ts,CbF,M01,M02,M1,M2,M001,F211,m01m1,m01m2,m02m1,m02m2,outp)

    return outp

end 
