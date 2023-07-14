# greens_optimization_comparison
Comparison of different ODE, numerical differentiation, and numerical integration methods for 2d integration with Green's theorem

This code compares different methods for approximating the line integral required for the Green's theorem approach for optimization.
We consider Euler's method and RK4 for approximating the integrand.
We consider 1st, 2nd, and 4th order FDM for approximating the derivative of the parametrized boundary.
We consider left hand rule, trapezoidal rule, and Gauss quadrature for the quadrature of the line integral.

We consider two examples with the integrands f(x,y)=x^2+y^2 and f(x,y)=cos(1/2*( x^2 + y^2)) on the "C" shaped domain.

Run the first example with:

include("greens_comp.jl")

greens_comp()

Run the second example with:

include("greens_comp2.jl")

greens_comp()
