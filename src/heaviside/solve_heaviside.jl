# Load solution procedures
include("functions_heaviside.jl")

"""
   Solve the sytem of second order ODEs 

   ``M A(t) + C V(t) + K U(t) = F(t)``

   with 

   U(t0) = U0

   V(t0) = V0

   for second order polynomials multiplied by Heaviside functions
    
   ``f_j = (c_jk0 + c_jk1 t + c_jk2 t^2)*H(t-t_jk)``.  

   Loading is informed by using a dictionary  ``load_data = Dict{Int64,Vector{Float64}}()``

   ``load_data[j] = [c_j00; c_j01; c_j02; t_j0; ... ; c_j(nk)0; c_j(nk)1; c_j(nk)2; t_j(nk))]``

   were nk is the number of coefficients used to represent the loading at DOF j.
   
   The output is 

   y(t)  => complete solution


"""
function Solve_heaviside(M::AbstractMatrix{T}, C::AbstractMatrix{T},K::AbstractMatrix{T},
                         U0::AbstractVector{T},V0::AbstractVector{T}, 
                         load_data::OrderedDict{Int64,Vector{Float64}}; tspan=(0.0,10.0),t0=0.0) where T

    # Basic assertions
    @assert tspan[1] < tspan[2] "Solve_heaviside:: Initial time must be smaller than the final time"
    @assert tspan[1] <= t0 <= tspan[2] "Solve_heaviside:: t0 must be in tspan"
    @assert t0==0 "Solve_heaviside:: t0 must be 0.0 by now"

    # Evaluate F211 
    chol = cholesky(M)
    Kb = chol\K
    Cb = chol\C
    F211 = 0.5*Cb + 0.5*sqrt(Cb^2 - 4*Kb)

    # Pre-evaluate matrices needed to compute the permanente solution
    CbF = Cb - F211
    M01 = CbF^(-1)
    M02 = CbF^(-2)
    M03 = CbF^(-3)
    M1 = F211^(-1)
    M2 = F211^(-2)
    M3 = F211^(-3)
    M001 = (CbF)^(-1)

    # Pre-process 
    sol_j = Process_heaviside(M,load_data)

    # Permanent solution for a given time
    yp(t) = y_permanent_heaviside(t,sol_j,load_data,CbF, M01, M02, M03, M1, M2, M3, M001, F211)

    # Complete response
    y(t) = yp(t)

    # Return the complete solution
    return y

end
