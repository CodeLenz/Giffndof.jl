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
                     load_data::OrderedDict{Int64,Vector{Float64}}; t0=0.0) where T

   
    # Basic assertions
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

    # Permanent solution for a given time
    yp(t) = y_permanent_dirac(t,sol_j,load_data,F211,CbF)

    # Complete response
    y(t) = yp(t)

    # Return complete response
    return y

end
