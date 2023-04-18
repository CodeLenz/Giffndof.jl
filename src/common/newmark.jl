"""
Solve the transient problem M(x)A(x,t) + C(x)V(x,t) + K(x,t)U(x,t) = F(t), 
using Newmark-beta method. 

    Solve_newmark(M::AbstractMatrix,C::AbstractMatrix,K::AbstractMatrix, f!::Function, gls::Matrix{Int64}, 
                  ts::Tuple{Float64, Float64}, Δt::Float64,
                  verbose=false;
                  U0=Float64[], V0=Float64[],
                  β=1/4, γ=1/2)

where 

    M, C and K are given n x n matrices

    ts is Tupple with initial and end time (Ti,Tf)  

    Δt is (fixed) time step  

    verbose is false or true

    U0 and V0 are the initial conditions  

    β and γ are the parameters of the Newmark method
 
    f!(t,F) must be a function of t  and F where F is n x 1,
                Example: 
    
                function f!(t,F)
                            P = [0.0;10.0]
                            F.= cos(2*t)*P 
                end   

Return three arrays of size n x nt, where  nt is the number of time steps (length of t0:Δt:tf)

    A_U displacements

    A_V velocities

    A_A accelerations
    
    A_t is a vector of size nt x 1 with discrete times

"""
function Solve_newmark(M::AbstractMatrix,C::AbstractMatrix,K::AbstractMatrix, f!::Function, 
                      ts::Tuple{Float64, Float64}, Δt::Float64;
                      U0=Float64[], V0=Float64[], β=1/4, γ=1/2)


    #
    #                              Initialization
    #             

    # Tspan    
    t0,tf = ts
    tspan = t0:Δt:tf-Δt

    # Number of time steps
    nt = length(tspan)

    # Problem's size
    nfull = size(M,1)

    # Newmark operator
    MN = M .+ β*K*Δt^2 .+ γ*C*Δt

    # Check initial conditons
    if isempty(U0)
        resize!(U0,nfull)
        fill!(U0,0.0)
    else
        length(U0)==nfull || throw("Newmark::U0 should be $nfull")
    end

    if isempty(V0)
        resize!(V0,nfull)
        fill!(V0,0.0)
    else
        length(V0)==nfull || throw("Newmark::V0 should be $nfull")
    end

    # Force vector
    F = zeros(nfull)

    # Initial acceleration
    f!(0.0,F)
    A0 = M\(F- K*U0 -C*V0)

    # Arrays to monitor the solution. The number of time points is 
    # tspan / dt + 1.
    ndofs = nfull
    A_t = Vector{Float64}(undef,nt+1)
    A_U = Array{Float64}(undef,nt+1,ndofs)
    A_V = Array{Float64}(undef,nt+1,ndofs)
    A_A = Array{Float64}(undef,nt+1,ndofs)

    # Store initial values (time zero)
    A_t[1]    = t0
    A_U[1,:] .= U0[:]
    A_V[1,:] .= V0[:]
    A_A[1,:] .= A0[:]

    # Main loop
    count = 2
    for t in tspan

            f!(t+Δt,F)
            b = F .- K*U0 .-(C .+Δt*K)*V0 .- (C*Δt*(1-γ) .+ K*(1/2-β)*Δt^2)*A0
            A = MN\b
            V = V0 .+ Δt*( (1-γ)*A0 .+ γ*A )
            U = U0 .+ Δt*V0 .+ ( (1/2-β)*A0 .+ β*A )*Δt^2
            
           # Store values at t+Δt
            A_t[count]    = t + Δt
            A_U[count,:] .= U[:]
            A_V[count,:] .= V[:]
            A_A[count,:] .= A[:]
            count += 1

            U0 .= U
            V0 .= V
            A0 .= A
    end

    # Return the values 
    return A_U, A_V, A_A, A_t

end


