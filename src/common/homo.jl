
"""
 Evaluate constants C1 and C2 due to intial conditions

 t0   -> time to impose the initial conditions

 F211 -> Constant matrix

 FCb  -> Constant matrix F211 - Cb

 Cb2F -> Constant matrix Cb - 2F211

 U0   -> Initial conditions - user provided real valued vector #

 V0   -> Initial conditions - user provided real valued vector 

 yp   -> Function returning permament solution at t

 dyp  -> Function returning the derivative of the permament solution at t

 return two (complex valued) vectors C1 and C2

"""
function Evaluate_Cs(t0::Float64,F211::AbstractMatrix{T1},FCb::AbstractMatrix{T1},
                     Cb2F::AbstractMatrix{T1},
                     U0::Vector{T2}, V0::Vector{T2},
                     yp::Function, dyp::Function) where {T1,T2}

    if t0==0.0
       #
       # Specific case when t0 = 0.0 (Equation A.5)
       #
       C1, C2 = Evaluate_Cs_t00(F211,FCb,Cb2F,U0,V0,yp(t0),dyp(t0))
    else
        #
        # General case for t0 different than 0.0 (Equations A.8 and A.9)
        #
       C1, C2 = Evaluate_Cs_general(t0,F211,FCb,U0,V0,yp(t0),dyp(t0))
    end

    return C1, C2
end


"""
 Evaluate Constants for homogeneous solution - Equations A.8 and A.9

 PARTICULAR CASE WHEN t0 = 0s


 F211 -> Constant matrix

 FCb  -> Constant matrix F211 - Cb

 Cb2F -> Constant matrix Cb - 2F211

 U0   -> Initial conditions - user provided real valued vector 

 V0   -> Initial conditions - user provided real valued vector 

 yp   -> Function returning permament solution at t

 dyp  -> Function returning the derivative of the permament solution at t

 return two (complex valued) vectors C1 and C2

"""
function Evaluate_Cs_t00(F211::AbstractMatrix{T1},FCb::AbstractMatrix{T1},
                     Cb2F::AbstractMatrix{T1},
                     U0::Vector{T2}, V0::Vector{T2},
                     Y0::Vector{T3}, dY0::Vector{T3}) where {T1,T2,T3}

    C2 = (Cb2F)\(V0 .-dY0 .- FCb*(U0.-Y0))
    C1 = U0 .- Y0 .- C2

    return C1, C2

end

"""
 Evaluate Constants for homogeneous solution - Equation A.5

 GENERAL CASE


 t0   -> time to impose the initial conditions

 F211 -> Constant matrix

 FCb  -> Constant matrix F211 - Cb

 Cb2F -> Constant matrix Cb - 2F211

 U0   -> Initial conditions - user provided real valued vector 

 V0   -> Initial conditions - user provided real valued vector 

 yp   -> Function returning permament solution at t

 dyp  -> Function returning the derivative of the permament solution at ts

 return two (complex valued) vectors C1 and C2
"""
function Evaluate_Cs_general(t0::Float64,F211::AbstractMatrix{T1},FCb::AbstractMatrix{T1},
                     U0::Vector{T2}, V0::Vector{T2},
                     Y0::Vector{T3}, dY0::Vector{T3}) where {T1,T2,T3}


    # Problem size
    n = length(U0)

    # Pre-compute exp(FCb*t0)
    expFCb = exp(FCb*t0)

    # Pre-compute exp(-F211*t0)
    PC3 = exp(-F211*t0)

    # Augmented linear system                 
    A = [  expFCb        PC3 ;
         FCb*expFCb   -F211*PC3 ]

    # Solve the augmented system
    Cs = A\[U0 .- Y0 ; V0 .- dY0]
    
    return Cs[1:n], Cs[n+1:end]

end

"""
 Evaluate homogeneous response at time t

 t      -> time

 F211   -> Constant matrix

 FCb    -> Constant matrix = F211 - Cb

 C1, C2 -> Constant vectors evaluated by Evaluate_Cs()


 return (complex valued) vector callable function y_h(t)
"""
function y_homo(t::Float64,F211::AbstractMatrix{T1}, FCb::AbstractMatrix{T2},
                C1::Vector{T3}, C2::Vector{T3}) where {T1,T2,T3}

    exp(-F211*t)*C2 + exp(FCb*t)*C1 

end


