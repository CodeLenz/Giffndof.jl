# Load solution procedures
include("functions_polynomial.jl")

"""
   Solve the sytem of second order ODEs 

   ``M A(t) + C V(t) + K U(t) = F(t)``

   with 

   U(t0) = U0

   V(t0) = V0

   for polynomial excitations  ``f_j = c_jk (t-tj)^k``.  

   Loading is informed by using a dictionary  ``load_data = Dict{Int64,Vector{Float64}}()``

   ``load_data[j] = [c_j1; c_j2, ... t_j]``

   were nk is the number of coefficients used to represent the 
   loading at DOF j.
   
   The output are three functions of t

   y(t)  => complete solution

   yh(t) => homogeneous solution

   yp(t) => permanent solution
     

"""
function Solve_polynomial(M::AbstractMatrix{T}, C::AbstractMatrix{T},K::AbstractMatrix{T},
                          U0::AbstractVector{T},V0::AbstractVector{T}, 
                          load_data::OrderedDict{Int64,Vector{Float64}}; t0=0.0) where T

   
    # Evaluate F211 - Equation 65
    chol = cholesky(Symmetric(M))
    Kb = chol\K
    Cb = chol\C
    F211 = 0.5*Cb + 0.5*sqrt(Cb*Cb - 4*Kb)

    # Evaluate FC (Auxiliary matrix to avoid repeated computation)
    FCb = F211 .- Cb

    # Evaluate Cb2F (Auxiliary matrix to avoid repeated computation)
    Cb2F = Cb .- 2*F211

    # Pre-process - cache some expensive operations
    sol_jklp = Process_polynomial(M, Kb, F211, load_data)

    # Permanent solution for a given time
    yp(t) = y_permanent_polynomial(t,sol_jklp,load_data)

    # Derivative of yp for a given time
    dyp(t) = dy_permanent_polynomial(t,sol_jklp,load_data)

    # Evaluate constants C1 and C2 - Appendix A
    C1, C2 = Evaluate_Cs(t0,F211,FCb,Cb2F,U0,V0,yp,dyp)

    # Homogeneous solution at t - Equation 49
    yh(t) = y_homo(t,F211,FCb,C1,C2)

    # Complete response
    y(t) = yp(t) + yh(t)

    # Return the complete,the homogeneous and the particular solutions
    return y, yh, yp

end
