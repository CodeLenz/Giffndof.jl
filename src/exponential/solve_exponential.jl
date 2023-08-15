
# Load solution procedures
include("functions_exponential.jl")

"""
   Solve the sytem of second order ODEs 

   ``M A(t) + C V(t) + K U(t) = F(t)``

   with 

   U(t0) = U0

   V(t0) = V0

   for exponential excitations  ``f_j = c_jk exp (im w_jk t + im phi_jk)``.  

   Loading is informed by using a dictionary  ``load_data = Dict{Int64,Vector{ComplexF64}}()``

   ``load_data[j] = [c_j1; w_j1; phi_j1; ....; c_jnk; w_jnk; phi_jnk]``

   were nk is the number of exponentials used to represent the 
   loading at DOF j.
   
   The output are three functions of t

   y(t)  => complete solution

   yh(t) => homogeneous solution

   yp(t) => permanent solution
     

"""
function Solve_exponential(M::AbstractMatrix{T}, C::AbstractMatrix{T},K::AbstractMatrix{T},
                           U0::AbstractVector{T},V0::AbstractVector{T}, 
                           load_data::Dict{Int64,Vector{ComplexF64}}; t0=0.0) where T

    
    # Pre-process - cache some expensive operations
    sol_jk, beta_jk, cphi_jk = Process_exponential(M, C, K, load_data)

    # Permanent solution for a given time t - Equation 120
    yp(t) = y_permanent_exponential(t,sol_jk,beta_jk,cphi_jk)

    # Derivative of yp(t) for a given time t - Equation 121
    dyp(t) = dy_permanent_exponential(t,sol_jk,beta_jk,cphi_jk)

    # Evaluate F211 - Equation 65
    chol = cholesky(Symmetric(M))
    Kb = chol\K
    Cb = chol\C
    F211 = 0.5*Cb + 0.5*sqrt(Cb^2 - 4*Kb)

    # Evaluate FC (Auxiliary matrix to avoid repeated computation)
    FCb = F211 .- Cb

    # Evaluate Cb2F (Auxiliary matrix to avoid repeated computation)
    Cb2F = Cb .- 2*F211

    # Evaluate constants C1 and C2 - Appendix A
    C1, C2 = Evaluate_Cs(t0,F211,FCb,Cb2F,U0,V0,yp,dyp)

    # Homogeneous solution at t - Equation 49
    yh(t) = y_homo(t,F211,FCb,C1,C2)

    # Complete response at t 
    y(t) = yp(t) + yh(t)

    # Return the complete,the homogeneous and the particular solutions
    return y, yh, yp

end
