# AutoTortoise.jl

`AutoTortoise.jl` is a small Julia package prototype for horizon-aware tortoise-coordinate construction for static spherically symmetric black-hole metrics whose metric function is of Laurent-polynomial form

```math
f(r)=\sum_n a_n r^n,
```

where the integer exponent `n` may be negative, zero, or positive.

The package is intended as infrastructure for QNM, CSM, and CLD calculations.  It automates the pipeline

```text
MetricFunction -> Horizons -> TortoiseMap -> InverseMap -> V(r(r_*))
```

for Laurent-polynomial metric functions.

## Current scope

Version `0.1.4` supports:

- Laurent-polynomial metric functions.
- Conversion to a rational integral via `P(r)=r^m f(r)`.
- Polynomial root finding with a companion matrix.
- Positive-real horizon detection.
- Simple-pole partial fractions for `1/f(r)`.
- Real tortoise maps on static intervals.
- Numerical inverse maps by bracketed bisection.
- Schwarzschild, Schwarzschild--dS, and Schwarzschild--AdS helper constructors.
- An experimental complex inverse map for CSM-style continuation.

Near-extremal or exactly degenerate horizons are detected but not yet treated as a fully reliable repeated-pole problem.

## Installation

```julia
] add AutoTortoise
```

## Basic usage

```julia
using AutoTortoise

f = schwarzschild_ds(M=1.0, Lambda=0.08)
h = horizons(f)
I = static_region(h; prefer=:finite)
tm = tortoise_map(f, I)

r = (I.left + I.right)/2
x = tortoise(tm, r)
r_back = inverse_tortoise(tm, x)
```

For a custom Laurent-polynomial metric:

```julia
f = MetricFunction(-1 => -2.0, 0 => 1.0, 2 => -0.08/3)
```

This corresponds to

```math
f(r)=1-\frac{2}{r}-\frac{0.08}{3}r^2.
```

## AdS-like finite boundary convention

For Schwarzschild--AdS,

```julia
f = schwarzschild_ads(M=1.0, L=10.0)
h = horizons(f)
I = static_region(h; prefer=:exterior)
tm = tortoise_map(f, I; convention=:boundary_zero)
```

The convention `:boundary_zero` normalizes the finite tortoise-coordinate boundary at infinity to `r_* = 0` for standard AdS-like rational metrics.

## Mathematical structure

For

```math
f(r)=\sum_n a_n r^n,
```

choose

```math
m=\max(0,-n_{\min}),
\qquad
P(r)=r^m f(r).
```

Then

```math
\frac{1}{f(r)}=\frac{r^m}{P(r)}.
```

For simple roots `r_i` of `P(r)`, the pole residues are

```math
A_i=\frac{r_i^m}{P'(r_i)}=\frac{1}{f'(r_i)}
```

when `r_i` is a nonzero simple horizon.  The tortoise coordinate is represented as

```math
r_*(r)=\sum_i A_i \log(r-r_i)+\int^r Q(\bar r)\,d\bar r+C,
```

where `Q` is the polynomial quotient produced by rational division.

On a real physical interval, real horizon terms are evaluated using the real branch `log(abs(r-r_i))`.

## Roadmap

Planned extensions:

- Fully reliable repeated-pole support for extremal and near-extremal horizons.
- BigFloat and ArbNumerics support for near-degenerate roots.
- Explicit branch tracking for complex-scaled rays.
- Grid generation and potential composition utilities for CSM calculations.
- Documentation examples for Reissner--Nordström, SdS, and AdS cases.
