
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
# Split  A*sin( f t + d)
#
# where f is in Hz and d in degrees
#
# into two exponentials
#
# A*i/2 * exp (-i(w + p)) 
#
# and
#
# -A*i/2 * exp (+i(w + p))
#
# where w is in rad/s and p in rad
#
function Split_sin(A,f,p)

    # Complex amplitudes
    c_j1 =  im*A/2
    c_j2 = -im*A/2

    # Frequencies
    w_j1 = -2*pi*f*1im
    w_j2 =  2*pi*f*1im

    # Phases
    p_j1 = -(p*(pi/180))*1im
    p_j2 =  (p*(pi/180))*1im

    # Return in the same sequence used in exponential
    [c_j1; w_j1; p_j1; c_j2; w_j2; p_j2]

end

#
# Split  A*cos( f t + d)
#
# where f is in Hz and d in degrees
#
# into two exponentials
#
# A*i/2 * exp (i(w + p)) 
#
# and
#
# A*i/2 * exp (-i(w + p))
#
# where w is in rad/s and p in rad
#
function Split_cos(A,f,p)

    # Complex amplitudes
    c_j1 =  A/2
    c_j2 = A/2

    # Frequencies
    w_j1 = 2*pi*f*1im
    w_j2 =  -2*pi*f*1im

    # Phases
    p_j1 = (p*(pi/180))*1im
    p_j2 =  -(p*(pi/180))*1im

    # Return in the same sequence used in exponential
    [c_j1; w_j1; p_j1; c_j2; w_j2; p_j2]

end