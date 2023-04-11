# Giffndof

[![Build Status](https://github.com/CodeLenz/Giffndof.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/CodeLenz/Giffndof.jl/actions/workflows/CI.yml?query=branch%3Amain)


# Generalized Integrating Factor for NDOFs

This repository contains the computer implementation of the solution procedures developed in <Link>

As this package is not yet registered, you must install it by using

```julia
using Pkg
Pkg.add("git@github.com:CodeLenz/Giffndof.jl.git")
```

Some examples of the manuscript are used to exemplify the package

## Example 1

```julia
#
# Solve a 3DOF problem subjected to an harmonic force 3sin(4t)
# in the second DOF.
#
function Example_exponential(;tspan = (0.0, 10.0), dt=0.01, t0 = 0.0)

    # Mass matrix
    M = [2.0 0.0 0.0 ;
         0.0 2.0 0.0 ;
         0.0 0.0 1.0 ]

    # Stiffness matrix
    K = [6.0 -4.0  0.0 ;
        -4.0  6.0 -2.0 ;
         0.0 -2.0  6.0]*1E2

    # Damping matrix
    C = 1E-2*K

    # Initial Conditions
    U0  = [0.0; 0.0; 0.0]
    V0  = [0.0; 0.0; 0.0]

    #----------------------------- g_2(t) = 3*sin( 4 t) -----------------------------#

    # Amplitude
    ampl = 3.0

    # Angular frequency
    ws = 4.0

    # Split the ampl*sin(ws t) into two exponentials

    # with apmplitudes
    c_21 =  ampl*im/2
    c_22 = -ampl*im/2

    # and angular frequencies
    w_21 = -ws
    w_22 =  ws

    # Create a dictionary. Each key corresponds to the DOF (j)
    # such that
    # load_data[j] = [c_j1; w_j1; ....; c_jnk; w_jnk]
    # were nk is the number of exponentials used to represent the 
    # loading at DOF j
    #
    load_data = Dict{Int64,Vector{ComplexF64}}()

    # For our example, DOF=2 and we have two exponentials
    load_data[2] = [c_21; w_21; c_22; w_22]

    # Main function -> solve the problem
    y, yh, yp = Solve_exponential(M,C,K,U0,V0,load_data,tspan=tspan,t0=t0)

    # Discrete times to make the plot
    tt = tspan[1]:dt:tspan[2]
      
    # Reshape to plot
    ndofs = size(M,1)
    yy = reshape([real(y(t))[k] for k=1:ndofs for t in tt],length(tt),ndofs)

    # Plot
    display(plot(tt,yy))

    # Return y, yh and yp
    return y, yh, yp

end

```

## Example 2

```julia
function Example_polynomial(;tspan = (0.0, 10.0), dt=0.01, t0 = 0.0)

    # Mass matrix
    M = [2.0 0.0 0.0 ;
         0.0 2.0 0.0 ;
         0.0 0.0 1.0 ]

    # Stiffness matrix
    K = [6.0 -4.0  0.0 ;
        -4.0  6.0 -2.0 ;
         0.0 -2.0  6.0]*1E2

    # Damping matrix
    C = 1E-2*K

    # Initial Conditions
    U0  = [0.0; 0.0; 0.0]
    V0  = [0.0; 0.0; 0.0]

    # Loading
    load_data = OrderedDict{Int64,Vector{Float64}}()

    # 10t - t^2 at DOF 2
    #        DOF     c20  c21   c22   t2
    load_data[2] = [0.0 ; 0.0; 10.0; -1.0]


    #  Main function -> solve the problem
    y, yh, yp = Solve_polynomial(M,C,K,U0,V0,load_data,tspan=tspan,t0=t0)

    # Discrete times to make the plot
    tt = tspan[1]:dt:tspan[2]
       
    # Reshape to plot
    ndofs = size(M,1)
    yy = reshape([real(y(t))[k] for k=1:ndofs for t in tt],length(tt),ndofs)
 
    # Plot
    display(plot(tt,yy))
 
    # Return y, yh and yp
    return y, yh, yp

end
```

## Example 3

```julia
function Example_dirac(;tspan = (0.0, 10.0), dt=0.01, t0 = 0.0)


    # Mass matrix
    M = [2.0 0.0 0.0 ;
         0.0 2.0 0.0 ;
         0.0 0.0 1.0 ]

    # Stiffness matrix
    K = [6.0 -4.0  0.0 ;
        -4.0  6.0 -2.0 ;
         0.0 -2.0  6.0]*1E2

    # Damping matrix
    C = 1E-2*K

    # Initial Conditions
    U0  = [0.0; 0.0; 0.0]
    V0  = [0.0; 0.0; 0.0]

    #
    # Loading 
    #
    # g_2(t) = delta(t-1) - delta(t-5)
    #
    load_data = OrderedDict{Int64,Vector{Float64}}()
    #
    #               c_20  t_20  c_21  t_21
    load_data[2] = [1.0 ; 1.0; -1.0 ; 5.0]

    #  Main function -> solve the problem
    y = Solve_dirac(M,C,K,U0,V0,load_data,tspan=tspan,t0=t0)

    # Discrete times to make the plot
    tt = tspan[1]:dt:tspan[2]
       
    # Reshape to plot
    ndofs = size(M,1)
    yy = reshape([real(y(t))[k] for k=1:ndofs for t in tt],length(tt),ndofs)
 
    # Plot
    display(plot(tt,yy))
 
    # Return y
    return y

end
```

## Example 4
```julia
function Example_heaviside(;tspan = (0.0, 10.0), dt=0.01, t0 = 0.0)

    # Mass matrix
    M = [2.0 0.0 0.0 ;
         0.0 2.0 0.0 ;
         0.0 0.0 1.0 ]

    # Stiffness matrix
    K = [6.0 -4.0  0.0 ;
        -4.0  6.0 -2.0 ;
         0.0 -2.0  6.0]*1E2

    # Damping matrix
    C = 1E-2*K

    # Initial Conditions
    U0  = [0.0; 0.0; 0.0]
    V0  = [0.0; 0.0; 0.0]

    #
    # Loading (1 + 0*t + 0*t^2) H(t-1) - (1 + 0*t + 0*t^2)H(t-5)
    #
    load_data = OrderedDict{Int64,Vector{Float64}}()

    #   c_j00 c_j01 c_j02  t_jk .... c_j(nk)0 c_j(nk)1 c_j(nk)2 t_j(nk)
    load_data[2] = [1.0; 0.0; 0.0; 1.0 ; -1.0; 0.0; 0.0; 5.0 ]

    #  Main function -> solve the problem
    y = Solve_heaviside(M,C,K,U0,V0,load_data,tspan=tspan,t0=t0)

    # Discrete times to make the plot
    tt = tspan[1]:dt:tspan[2]
        
    # Reshape to plot
    ndofs = size(M,1)
    yy = reshape([real(y(t))[k] for k=1:ndofs for t in tt],length(tt),ndofs)
  
    # Plot
    display(plot(tt,yy))
  
    # Return y
    return y 

end
```