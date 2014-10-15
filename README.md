# FiniteDifferenceDerivatives

[![Build Status](https://travis-ci.org/pwl/FiniteDifferenceDerivatives.jl.svg?branch=master)](https://travis-ci.org/pwl/FiniteDifferenceDerivatives.jl)

Usage example
=============

Function `fdd(k,n,x,f)` computes the `k`-th derivative of `f` on a
nonuniformly spaced mesh `x` using `n`-point stencli.

```
using FiniteDifferenceDerivatives
npts = 11
x = linspace(0,1,npts)
f = x.^2
df = fdd(1,3,x,f) # Compute the first derivative using three point stencil
```
