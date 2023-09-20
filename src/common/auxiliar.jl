
#
# Avoid problems with Sparse Matrices in exp(A)
#
import LinearAlgebra.exp
function exp(A::SparseMatrixCSC)
    exp(Array(A))
end

import LinearAlgebra.sqrt
function sqrt(A::SparseMatrixCSC)
    sqrt(Array(A))
end


#
# Split  A*sin(wt + p) into two exponentials
#
#
# A*i/2 * exp (-i(w + p)) 
#
# and
#
# -A*i/2 * exp (+i(w + p))
#
#
function Split_sin(A,w,p)

    # Complex amplitudes
    c_j1 =  im*A/2
    c_j2 = -im*A/2

    # Frequencies
    w_j1 = -w*1im
    w_j2 =  w*1im

    # Phases
    p_j1 = -p*1im
    p_j2 =  p*1im

    # Return in the same sequence used in exponential
    [c_j1; w_j1; p_j1; c_j2; w_j2; p_j2]

end

#
# Split  A*cos(wt + p) into two exponentials
#
#
# A*i/2 * exp (i(w + p)) 
#
# and
#
# A*i/2 * exp (-i(w + p))
#
#
function Split_cos(A,w,p)

    # Complex amplitudes
    c_j1 = A/2
    c_j2 = A/2

    # Frequencies
    w_j1 = w*1im
    w_j2 =  -w*1im

    # Phases
    p_j1 = p*1im
    p_j2 =  -p*1im

    # Return in the same sequence used in exponential
    [c_j1; w_j1; p_j1; c_j2; w_j2; p_j2]

end