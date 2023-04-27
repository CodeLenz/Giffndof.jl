"""
 Evaluate the coefficients a used to compute coefficients c.

 Continuous version (load is a function of time)

 Inputs:

 g is the reference (known) function 

 tl is the initial time at the interval

 tL is the final time at the interval

 Outputs:

 a0 coefficient 

"""
function Evaluate_coefs_aH0(g::Function, tl::Float64, tL::Float64)

    # Testing purposes
    tL > tl || error("Evaluate_coefs_aH0:: tl must be smaller than tL")

    # Values of g(t) at the extremes of the interval
    gl = g(tl)
    gL = g(tL)

    # Evaluate dt for this interval
    dt = tL-tl

    # Integral of g in the interval. 
    # Evalute using quadgk package.
    integral = quadgk(g,tl,tL)[1] 

    # Coeficient a0
    a0 = integral/dt 

    # Return a0
    return a0

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

"""
function Evaluate_coefs_aH0(interval::Int64,g::Vector{T}, Ts::T1) where {T,T1}

    # Recover the times
    tl = Ts[interval]
    tL = Ts[interval+1]

    # Values of g at the extremes of the interval
    gl = g[interval]
    gL = g[interval+1]

    # Evaluate dt for this interval
    dt = tL-tl

    # Integral of g in the interval. 
    a0 = (gl + 0.5*(gL-gl)) 

    # Return a0
    return a0

end


"""
 Evaluate the coefficients c. Continuous version.

 Inputs:

 g is the reference (known) function 

 Ts a list or StepRange with discrete times
 
 Outputs:

 cj0 coefficient 
 
"""
function Evaluate_coefs_cH0(g::Function, Ts::T) where T

   # Number of intervals 
   n = length(Ts)-1

   # Allocate both output vectors
   cj0 = zeros(n)

   # Intial time
   tl = Ts[1]

   # Loop over the intervals
   for interval=1:n

       # Final time at this interval
       tL = Ts[interval+1]
       
       # Evaluate coef a0 for this interval
       a0 = Evaluate_coefs_aH0(g,tl,tL)

       # Evaluate coef c
       cj0[interval] = a0 - sum(cj0[1:interval-1])
    
       # Update tl
       tl = tL

   end # interval

   # Return coefs cj0
   cj0

end

"""
 Evaluate the coefficients c. Discrete

 Inputs:

 g is a vector with discrete values

 Ts is a vector of discrite times or a StepRange

 Outputs:

 cj0 coefficient 
 
"""
function Evaluate_coefs_cH0(g::Vector{T}, Ts::T1) where {T,T1}

    # Basic asseertion
    length(g)==length(Ts) || error("Evaluate_coefs_cH0:: vectors g and Ts must have the same length")

   # Number of intervals 
   n = length(Ts)-1

   # Allocate output vector
   cj0 = zeros(n)

   # Loop over the intervals
   for interval=1:n
 
       # Evaluate coef a0 for this interval
       a0 = Evaluate_coefs_aH0(interval,g,Ts)

       # Evaluate coefs c
       cj0[interval] = a0 - sum(cj0[1:interval-1])
     
   end # interval

   # Return the coefficients
   cj0

end


#
# Evaluate all coefficients cj0 and cj1 
# and store in a OrderedDict
#
function Generate_Dict_cH0(load_data::OrderedDict{Int64,T},Ts::T1) where {T,T1}

    # Create the dictionary with c's for every DOF j informed in load_data
    dict_c = OrderedDict{Int64,Vector{Float64}}()

    # Loop over load_data
    for j in keys(load_data)

        # Evaluate the Coefficients c
        cj0 = Evaluate_coefs_cH0(load_data[j], Ts)

        # Store in the new dictionary
        dict_c[j] = cj0
        
    end #j

    # Return the dictionary with the c's for every DOF
    return dict_c

end


#
# Build the approximation for g(t) at t
#
"""
 Evaluate the approximation for a given function g(t) at t. 
 This approximation is built by using the Heaviside Series with cte values (order zero)
 at each interval. The purpose of this subroutine is for verification purposes only,
 since the evaluation of the response y⁽⁰⁾(t) is performed by using the coefficients
 c only.

 Inputs:

 t is the time to evaluate the approximation

 cj0 is a vector of coefficients (obtained by using Evaluate_coefs_c)

 Ts is a vector or StepRange with the discrete times used to evaluate cj0 and cj1

 Output:

 the approximation at
 
"""
function Evaluate_gtildeH0(t::Float64, cj0::Vector{Float64}, Ts::T) where T

   # Number of intervals in Ts
   n = length(Ts)-1

   # Basic test
   n>1 || error("Evaluate_gtildeH0:: you must inform at least two values in Ts")

   # Inicialize the output value
   saida = 0.0

   # We do not need to evaluate all intervals. Build a logical mask
   # where t is smaller than Ts (0)
   mask = t .<= Ts

   # Thus, the number of intervals to use in the loop is
   np = max(1,min(n+1 - sum(mask),n))

   # Loop over the valid intervals
   for interval=1:np

        # Initial time at the interval
        tl = Ts[interval]
          
        # Add the contribution for this interval
        saida += Heaviside(t,tl)*(cj0[interval])
   
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

M1 = F211^(-1)

M001 = Cb2F^(-1)

F211  -> constant matrix

m01m1 = M01*M1

outp         -> output vector (modified in place) 
"""
function y_permanent_HS0!(t::Float64,sol_j::AbstractMatrix,dict_c::OrderedDict{Int64,Matrix{Float64}},Ts::T,
                      CbF::AbstractMatrix,
                      M01::AbstractMatrix,  M1::AbstractMatrix,
                      M001::AbstractMatrix,
                      F211::AbstractMatrix, m01m1::AbstractMatrix, 
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
        cj0 = @view dict_c[j]
        
        # Loop over k just for valid intervals in Ts
        for k=1:np

            # Polynomial Coefficients
            c_jk0 = cj0[k]
            
            # offset (initial time of this interval)
            t_jk = Ts[k]
            
            # Common term
            t1 = t_jk - t

            # T5
            T5 =  c_jk0*m01m1

            # T6
            T6A = -(c_jk0)*M1
            T6 = exp(F211*t1)*M001*(T6A)

            # T7
            T7A = m01m1*(-c_jk0)   
            T7D = -M1*(c_jk0)

            # Final T7
            T7 = exp(CbF*t1)*(T7A  .+ T7D)

            # Add and use the cache
            outp .= outp .+ (T5 .+ T6 .+ T7)*sol_j[:,count]

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

M1 = F211^(-1)

M001 = Cb2F^(-1)

F211  -> constant matrix

m01m1 = M01*M1

outp         -> output vector (modified in place) 
"""
function y_permanent_HS0(t::Float64,sol_j::AbstractMatrix,dict_c::OrderedDict{Int64,Matrix{Float64}},Ts::T,
                      CbF::AbstractMatrix,
                      M01::AbstractMatrix, M1::AbstractMatrix,
                      M001::AbstractMatrix,
                      F211::AbstractMatrix, m01m1::AbstractMatrix)  where {T}



    # Number of DOFs
    ngls = size(sol_j,1)

    # Alocate the output vector
    outp = Vector{ComplexF64}(undef,ngls)

    # Call the driver
    y_permanent_HS0!(t,sol_j,dict_c,Ts,CbF,M01,M1,M001,F211,m01m1,outp)

    return outp

end 
