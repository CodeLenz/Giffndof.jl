#
# Solution for a Dirac (unitary impulse) excitation, Eq 150.
#
#  ``g_j(t) = sum_k c_jk Delta(t - t_jk)``
#
# where j is the DOF and k=1..nk is the number of base terms
#

"""
Pre-process data to avoid some computations when evaluating the 
permanent solution yp(t). Those computations are used many times,
so we evaluate them and store into a cache.

Inputs:

 M -> Constant matrix

 Kb -> Constant matrix ``Kb = M\\K``

 F211 -> Constant matrix

 ``load_data`` -> Dictionary with key j (gl) and data ``[c_j0 t_0 ...  c_jnk t_nk]``
 
Outputs:

 ``sol_jlp`` -> Matrix solutions of ``[C-2MF_211]\\e_j`` in each collumn

"""
function Process_dirac(M::AbstractMatrix{T}, C::AbstractMatrix{T},
                       F::AbstractMatrix{T1}, 
                       load_data::OrderedDict{Int64,Vector{Float64}}) where {T,T1}


    # Number of DOFs
    ngls = size(M,1)

    # ####################################################################
    #       Important fact [C-2MF]\e_j only depends on DOF j
    ######################################################################

    # Evaluate C-2MF 
    Coef = C .- 2*M*F

    # Lu decomposition
    LU = lu(Coef)

    # Number of entries to be stored in the cache
    ncol = length(load_data)

    # Cache 
    sol_j = zeros(ComplexF64,ngls,ncol)

    # Pre allocate e_j
    e_j = zeros(ngls)

    # Evaluate the cache
    count = 1
    for j in keys(load_data)

        # Data
        data = load_data[j]

        # Set position j to 1.0
        e_j[j] = 1.0

        # Solve [C-2MF]\e_j and store in the cache
        sol_j[:,count] .= LU\e_j

        # update count
        count += 1

        # Unset e_j
        e_j[j] = 0.0

    end #j

    # Return sol_j
    return sol_j

end

"""
Driver to evaluate the permanent solution at t 

Inputs:

t           -> time

``sol_j``  -> cache computed by Process_dirac

``load_data`` -> OrderedDict with loading information

F211  -> constant matrix

CbF -> constant matrix Cb - F211

outp         -> output vector (modified in place) 
"""
function y_permanent_dirac!(t::Float64,sol_j::AbstractMatrix{T},
                            load_data::OrderedDict{Int64,Vector{Float64}},
                            F211::AbstractMatrix{T},CbF::AbstractMatrix{T}, outp::Vector{T}) where T


    # Make sure outp is zero
    fill!(outp,zero(T))

    # Number of columns in the cache
    ncol = size(sol_j,2)

    #
    # Compute Eq. 179
    #
    for j in keys(load_data)

        # Data
        data = load_data[j]

        # Number of terms
        nk = Int(length(data)/2)

        # Counter (for j)
        count = 1

        # Loop over k
        for k=0:nk-1 

            # Coefficient
            c_jk = data[2*k+1]  

            # offset
            t_jk = data[2*k+2]  
            
            # Heaviside
            H_jk = Heaviside(t,t_jk)

            # Just evaluate if H is not null
            if H_jk > 0.0

                # Common term
                t1 = t_jk - t

                # First matrix
                M1 = exp(F211*t1)

                # Second matrix
                M2 = exp(CbF*t1)

                # Use the cache and add to the solution
                outp .= outp .+ c_jk*H_jk*(M1 .- M2)*sol_j[:,count]

            end # if H_jk > 0.0
            
        end #k

        # Update count
        count += 1

    end #j

end

"""
Evaluate the permanent solution at t 

Inputs:

t           -> time

``sol_j``  -> cache computed by Process_dirac

``load_data`` -> OrderedDict with loading information

F211  -> constant matrix

CbF -> constant matrix Cb - F211

Outputs:

outp         -> output vector
"""
function y_permanent_dirac(t::Float64,sol_j::AbstractMatrix{T},load_data::OrderedDict{Int64,Vector{Float64}},
                           F211::AbstractMatrix{T},CbF::AbstractMatrix{T}) where T

    # Number of DOFs
    ngls = size(sol_j,1)

    # Allocate output vector
    outp = Vector{T}(undef,ngls)

    # Call the driver
    y_permanent_dirac!(t,sol_j,load_data,F211,CbF,outp)

    # return
    return outp

end


