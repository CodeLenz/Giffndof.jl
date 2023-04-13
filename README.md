# Giffndof
![logo](Logo-Giff.png)


[![Build Status](https://github.com/CodeLenz/Giffndof.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/CodeLenz/Giffndof.jl/actions/workflows/CI.yml?query=branch%3Amain)


# Generalized Integrating Factor for NDOFs

This repository contains the computer implementation of the solution procedures developed in **to be included** 
for solving coupled systems of second order ODEs with constant coefficients

$$ M A(t) + C V(t) + K Y(t) = F(t) $$

with initial conditions

$$ Y(t_0) = U0 $$

and

$$ V(t_0) = V0 $$

where $t$ is the independent variable, $Y(t)$ a $n \times 1$ vector (dependent variables), $V$ its first   derivative with respect to $t$ and $A$ its second derivative. Matrices  $M$, $C$ and $K$ are $n \times n$. Vector $F(t)$ is informed by using a dictionary.
 
As this package is not yet registered, you must install it by using

```julia
using Pkg
Pkg.add("https://github.com/CodeLenz/Giffndof.jl.git")
```

***Disclaimer: This is a first version of this package and many of the optimizations discussed in the manuscript are not yet implemented***

There are four main methods being exported, depending on the type
of excitation:

```julia
y, yh, yp = Solve_exponential(M,C,K,U0,V0,load_data,t0=t0)
```

```julia
y, yh, yp = Solve_polynomial(M,C,K,U0,V0,load_data,t0=t0)
```

```julia
y, yh, yp = Solve_dirac(M,C,K,U0,V0,load_data,t0=t0)
```

```julia
y, yh, yp = Solve_heaviside(M,C,K,U0,V0,load_data,t0=t0)
```

where $y$ is the complete solution, $y_h$ the homogeneous solution and $y_p$ the permanent solution. Those solutions are functions and can be evaluated by simply passing a given time 

```julia
y(0.1) 
```

for example.


There is a specific way of informing non null entries in $F(t)$ for each type of excitation. 


## Exponentials

For forces described as a series of exponentials 

$$ f_j(t) = \sum_{k=1}^{n_k} c_{jk} \exp(i \omega_{jk} t + \phi_{jk}) $$

the user must inform the DOF $j$ as a key to a dictionary with entries given by (possible complex values) of $c_{jk}$ and $\omega_{jk}$

```julia
    load_data = Dict{Int64,Vector{ComplexF64}}()
```

Lets consider the first example in the reference manuscript

### Example 1

Consider a $3$ DOFs problem subjected to a force 

$$ f_2(t) = 3 \sin(4t) = 3\frac{i}{2}(\exp(-4it) - \exp(4it)) $$

such that the (complex) amplitudes are $c_{21}=3i/2$ and $c_{22}=-3i/2$ and the angular frequencies are $\omega_{21}=-4$ and $\omega_{22}=4$. Thus,

```julia
load_data[2] = [3*im/2; -3*im/2; -4.0; 4.0]
```

The complete example is 

```julia
using Giffndof
function Example_exponential(t0 =0.0)

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
    y, yh, yp = Solve_exponential(M,C,K,U0,V0,load_data,t0=t0)

    # Return the solutions for any t
    return y, yh, yp
    
 end   
```

One can generate the visualization for $y(t)$

```julia
  using Plots
  function Generate_plot(tspan = (0.0, 10.0), dt=0.01)

    # Call the example
    y, yh, yp = Example_exponential()

    # Discrete times to make the plot
    tt = tspan[1]:dt:tspan[2]
      
    # Reshape to plot
    ndofs = size(y(0.0),1)
    yy = reshape([real(y(t))[k] for k=1:ndofs for t in tt],length(tt),ndofs)

    # Plot
    display(plot(tt,yy))

end

```

## Polynomials

For forces described as a polynomial

$$ f_j(t) = \sum_{k=0}^{n_k} c_{jk} (t-t_j)^k $$

the user must inform the DOF $j$ as a key to a dictionary with entries given by of $c_{jk}$ and $t_j$

```julia
    load_data = Dict{Int64,Vector{Float64}}()
```

### Example 2

Consider a $3$ DOFs problem subjected to a force 

$$ f_2(t) = 10 t - t^2 $$

such that $t_2=0$,  $c_{20}=0$,  $c_{21}=10$,  $c_{22}=-1$. Thus

```julia
load_data[2] = [0.0; 0.0; 10.0; -1.0]
```

The complete example is 

```julia
using Giffndof
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
    #        DOF     t2   c20  c21    c22  
    load_data[2] = [0.0 ; 0.0; 10.0; -1.0]

    #  Main function -> solve the problem
    y, yh, yp = Solve_polynomial(M,C,K,U0,V0,load_data,t0=t0)

    # Return the solution
    return y, yh, yp

end

```

One can generate the visualization for $y(t)$

```julia
  using Plots  
  function Generate_plot(tspan = (0.0, 10.0), dt=0.01)

    # Call the example
    y, yh, yp = Example_polynomial()

    # Discrete times to make the plot
    tt = tspan[1]:dt:tspan[2]
      
    # Reshape to plot
    ndofs = size(y(0.0),1)
    yy = reshape([real(y(t))[k] for k=1:ndofs for t in tt],length(tt),ndofs)

    # Plot
    display(plot(tt,yy))

end
```

## Unitary impulse (Dirac's delta)

For forces described as a series of unitary impulses

$$ f_j(t) = \sum_{k=0}^{n_k} c_{jk} \delta(t-t_{jk}) $$

the user must inform the DOF $j$ as a key to a dictionary with entries given by of $c_{jk}$ and $t_{jk}$

```julia
    load_data = Dict{Int64,Vector{Float64}}()
```

### Example 3

Consider a $3$ DOFs problem subjected to two oposite unitary impulses at $t=1$ and $t=5$ s

$$ f_2(t) = \delta(t-1) - \delta(t-5) $$

such that $c_{20}=1.0$, $t_{20}=1$, $c_{21}=-1$ and $t_{21}=5.0$

```julia
load_data[2] = [1.0; 1.0; -1.0; 5.0]
```

The complete example is 

```julia
using Giffndof
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
    y, yh, yp = Solve_dirac(M,C,K,U0,V0,load_data,t0=t0)

    # Return the solution
    return y, yh, yp

 end
```
One can generate the visualization for $y(t)$

```julia
  using Plots
  function Generate_plot(tspan = (0.0, 10.0), dt=0.01)

    # Call the example
    y, yh, yp = Example_dirac()

    # Discrete times to make the plot
    tt = tspan[1]:dt:tspan[2]
      
    # Reshape to plot
    ndofs = size(y(0.0),1)
    yy = reshape([real(y(t))[k] for k=1:ndofs for t in tt],length(tt),ndofs)

    # Plot
    display(plot(tt,yy))

end
```

## Second order polynomials multiplied by Heavisides

For forces described as second order polynomials times heavisides

$$ f_j(t) = \sum_{k=0}^{n_k} (c_{jk0} + c_{jk1} t + c_{jk2} t^2) H(t-t_{jk}) $$

the user must inform the DOF $j$ as a key to a dictionary with entries given by of $c_{jk*}$ and $t_{jk}$

```julia
    load_data = Dict{Int64,Vector{Float64}}()
```

### Example 4

Consider a $3$ DOFs problem subjected to two oposite unitary steps at $t=1$ and $t=5$ s

$$ f_2(t) = H(t-1) - H(t-5) $$

such that $c_{200}=1$, $c_{201}=0$, $c_{202}=0$, $t_{20}=1$, $c_{210}=-1$, $c_{211}=0$, $c_{212}=0$, $t_{21}=5$

```julia
load_data[2] = [1.0; 0.0; 0.0; 1.0; -1.0; 0.0; 0.0; 5.0]
```

The complete example is 

```julia
using Giffndof
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
    y, yh, yp = Solve_heaviside(M,C,K,U0,V0,load_data,t0=t0)

    # Return the solution
    return y, yh, yp
    
 end
``` 

 One can generate the visualization for $y(t)$

```julia
  using Plots  
  function Generate_plot(tspan = (0.0, 10.0), dt=0.01)

    # Call the example
    y, yh, yp = Example_heaviside()

    # Discrete times to make the plot
    tt = tspan[1]:dt:tspan[2]
      
    # Reshape to plot
    ndofs = size(y(0.0),1)
    yy = reshape([real(y(t))[k] for k=1:ndofs for t in tt],length(tt),ndofs)

    # Plot
    display(plot(tt,yy))

end
```
