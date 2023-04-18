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
    w_j1 = -w
    w_j2 =  w

    # Phases
    p_j1 = -p
    p_j2 =  p

    # Return in the same sequence used in exponential
    [c_j1; w_j1; p_j1, c_j2; w_j2; p_j2]

end