"""
Pre-process data to avoid some computations when evaluating the 
permanent solution yp(t). Those computations are used many times,
so we evaluate them and store into a cache.

Inputs:

 M -> Constant matrix

 Kb -> Constant matrix ``Kb = M\\K``

 F211 -> Constant matrix

 ``load_data`` -> Dictionary with key j (gl) and data 
            ``[c_jk0 c_jk1 c_jk2 t_j1 ... c_j(nk)0 c_j(nk)1 c_j(nk)2 t_j(nk)]``
or 
            ``function`` 

    Data stored in the dictionary is not used in this subroutine.            

Outputs:

 ``sol_j`` -> Matrix solutions of ``M\\e_j`` in each collumn

"""
function Process_heaviside(M::AbstractMatrix{T},  
                           load_data::OrderedDict{Int64,T1}) where {T,T1}


    # Number of DOFs
    ngls = size(M,1)

    # ####################################################################
    #                Important [M]\e_j depends only on j
    ######################################################################

    # Cholesky decomposition
    chol = cholesky(M)

    # Number of cached solutions
    ncol = length(load_data)

    # Alocate the cache (Real valued matrix)
    sol_j = zeros(ngls,ncol)

    # Pre alocate e_j
    e_j = zeros(ngls)

    # Compute M\e_j
    count = 1
    for j in keys(load_data)

        # Set position j to 1.0
        e_j[j] = 1.0

        # Solve M\e_j
        sol_j[:,count] .= chol\e_j

        # update count
        count += 1

        # Unset e_j
        e_j[j] = 0.0

    end #j

    # Return sol_j
    return sol_j

end
