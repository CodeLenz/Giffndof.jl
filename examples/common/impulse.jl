
# Emulate dirac
function Impulse(t,t0,eps,A)

    a = 1/(2*eps)
    val = 0.0

    if t0-eps <= t <= t0+eps
        val = A*a*(1+cos(pi*(t-t0)/eps))
    end

    return val

end
