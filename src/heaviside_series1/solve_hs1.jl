# Load solution procedures
include("functions_hs1.jl")

"""
   Solve the sytem of second order ODEs 

   ``M A(t) + C V(t) + K U(t) = F(t)``

   with 

   U(t0) = U0

   V(t0) = V0

   when excitations are described by using first order Heaviside series built by 
   using known reference functions ``g(t)``.  

   Loading is informed by using a dictionary  ``load_data = Dict{Int64,Function}()``

   ``load_data[j] = [g(t)]``

   where ``g(t)`` is a function of time only.

   Ts is a vector of discrete times or a StepRange

   The output are three functions of t

   y(t)  => complete solution

   yh(t) => homogeneous solution

   yp(t) => permanent solution
     

"""
function Solve_HS1(M::AbstractMatrix{T}, C::AbstractMatrix{T},K::AbstractMatrix{T},
                   U0::AbstractVector{T},V0::AbstractVector{T}, Ts::T0,
                   load_data::OrderedDict{Int64,Function}; t0=0.0) where {T0,T}



    # Evaluate F211 
    chol = cholesky(M)
    Kb = chol\K
    Cb = chol\C
    F211 = 0.5*Cb + 0.5*sqrt(Cb^2 - 4*Kb)

    # Pre-evaluate matrices needed to compute the permanente solution
    CbF = Cb .- F211
    FCb = -CbF
    M01 = CbF^(-1)
    M02 = CbF^(-2)
    M1 = F211^(-1)
    M2 = F211^(-2)
    Cb2F = Cb .- 2*F211
    M001 = (Cb2F)^(-1)
    m01m1 = M01*M1
    m01m2 = M01*M2
    m02m1 = M02*M1
    m02m2 = M02*M2

    # Pre-process 
    sol_j = Process_heaviside(M,load_data)

    # Create dict_c
    dict_c = Generate_Dict_c(load_data,Ts)

    # Evaluate the permanent response
    yp(t) = y_permanent_HS1(t,sol_j,dict_c,Ts,CbF, M01, M02, M1, M2, M001, F211,  m01m1,
                            m01m2,m02m1,m02m2)
   
    # Evaluate constants C1 and C2 - Appendix A
    C1, C2 = Evaluate_Cs(t0,F211,FCb,Cb2F,U0,V0)

    # Homogeneous solution at t - Equation 49
    yh(t) = y_homo(t,F211,FCb,C1,C2)

    # Complete response at t 
    y(t) = yp(t) + yh(t)

    # Return the complete,the homogeneous and the particular solutions
    return y, yh, yp

end