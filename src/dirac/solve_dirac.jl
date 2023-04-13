# Load solution procedures
include("functions_dirac.jl")

"""
   Solve the sytem of second order ODEs 

   ``M A(t) + C V(t) + K U(t) = F(t)``

   with 

   U(t0) = U0

   V(t0) = V0

   for unitary impulses excitations  ``f_j = c_jk delta(t-t_jk)``.  

   Loading is informed by using a dictionary  ``load_data = Dict{Int64,Vector{Float64}}()``

   ``load_data[j] = [c_j1; t_j2; ... ; c_jnk;  t_jnk]``

   were nk is the number of coefficients used to represent the 
   loading at DOF j.
   
   The output is 

   y(t)  => complete solution


"""
function Solve_dirac(M::AbstractMatrix{T}, C::AbstractMatrix{T},K::AbstractMatrix{T},
                     U0::AbstractVector{T},V0::AbstractVector{T}, 
                     load_data::OrderedDict{Int64,Vector{Float64}}; tspan=(0.0,10.0),t0=0.0) where T

   
    # Basic assertions
    @assert tspan[1] < tspan[2] "Solve_dirac:: Initial time must be smaller than the final time"
    @assert tspan[1] <= t0 <= tspan[2] "Solve_dirac:: t0 must be in tspan"
    @assert t0==0 "Solve_dirac:: t0 must be 0.0 by now"

    # Evaluate F211 - Equation 65
    chol = cholesky(M)
    Kb = chol\K
    Cb = chol\C
    F211 = 0.5*Cb + 0.5*sqrt(Cb^2 - 4*Kb)
  
    # Evaluate CbF (Auxiliary matrix to avoid repeated computation)
    CbF =  Cb .- F211
  
    # Pre-process 
    sol_j = Process_dirac(M, C, F211, load_data)

    # Evaluate FC (Auxiliary matrix to avoid repeated computation)
    FCb = -CbF

    # Evaluate Cb2F (Auxiliary matrix to avoid repeated computation)
    Cb2F = Cb .- 2*F211

    # Evaluate constants C1 and C2 - Appendix A
    C1, C2 = Evaluate_Cs(t0,F211,FCb,Cb2F,U0,V0)

    # Homogeneous solution at t - Equation 49
    yh(t) = y_homo(t,F211,FCb,C1,C2)

    # Permanent solution for a given time
    yp(t) = y_permanent_dirac(t,sol_j,load_data,F211,CbF)

    # Complete response
    y(t) = yp(t) + yh(t)

    # Return complete response, homogeneous and particular
    return y, yh, yp

end