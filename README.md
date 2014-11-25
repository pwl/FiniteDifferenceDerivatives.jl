# FiniteDifferenceDerivatives

[![Build Status](https://travis-ci.org/pwl/FiniteDifferenceDerivatives.jl.svg?branch=master)](https://travis-ci.org/pwl/FiniteDifferenceDerivatives.jl)

Usage example
=============

Function `fdd!(df,f,x,k,n)` computes the `k`-th derivative of `f` on a
nonuniformly spaced mesh `x` using `n`-point stencli.  The result is
written to the array `df`.

```
using FiniteDifferenceDerivatives
x = linspace(0,1,11)
f = x.^2
df = zero(f)
fdd!(df,f,x,1,3) # Compute the first derivative using three point stencil
```

You can also use `fdd`, which allocates and returns the `df`
```
df = fdd(f,x,1,3)
```

It is also possible to compute the differentiation matrix with

```
D1 = fddmatrix(x,1,3)
using Base.Test
@test_approx_eq(D1*f,df)
```

To compute the second derivative at point `x0` given the function values `f` at
mesh points `x` you can call

```
x = linspace(0,1,11)
fddat(x.^2,x,2,0.5)
```

It is also possible to compute the derivative outside the interval
`[0,1]` but this might generate large errors for nonpolynomial
functions f

```
fddat(x.^2,x,1,2.0) # returns 4.000000002793968
fddat(1./(1.+x.^2),x,1,2.0) # returns 17.610649429727346
```

the true first derivative of `1/(1+x^2)` at `x=2.0` is `-0.16`.


TODO
====

- Add specialized functions for symmetric finite differences
