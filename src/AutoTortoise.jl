module AutoTortoise

using LinearAlgebra
using Logging

export MetricFunction,
       PhysicalInterval,
       Horizon,
       Horizons,
       PoleTerm,
       TortoiseMap,
       metric_value,
       metric_derivative,
       laurent_polynomial,
       polynomialized_metric,
       horizons,
       physical_intervals,
       static_region,
       tortoise_map,
       tortoise,
       tortoise_complex,
       inverse_tortoise,
       inverse_tortoise_complex,
       potential_at,
       schwarzschild,
       schwarzschild_lambda,
       schwarzschild_ds,
       schwarzschild_ads,
       surface_gravity

const DEFAULT_ROOT_TOL = 1.0e-8
const DEFAULT_IMAG_TOL = 1.0e-9

"""
    MetricFunction(coeffs)
    MetricFunction(pairs...)

Laurent-polynomial static metric function

```math
f(r) = sum_n a_n r^n,
```

where the exponents `n` are integers and may be negative, zero, or positive.

Examples
========

```julia
f = MetricFunction(Dict(-1 => -2.0, 0 => 1.0, 2 => -0.08/3))
f = MetricFunction(-1 => -2.0, 0 => 1.0, 2 => -0.08/3)
```
"""
struct MetricFunction
    coeffs::Dict{Int,ComplexF64}
    exponents::Vector{Int}
    nmin::Int
    nmax::Int
end

function MetricFunction(coeffs::AbstractDict{<:Integer,<:Number}; atol::Real=0.0)
    d = Dict{Int,ComplexF64}()
    for (n, a) in coeffs
        c = ComplexF64(a)
        if abs(c) > atol
            d[Int(n)] = c
        end
    end
    isempty(d) && error("MetricFunction requires at least one nonzero coefficient")
    exps = sort(collect(keys(d)))
    return MetricFunction(d, exps, first(exps), last(exps))
end

MetricFunction(pairs::Pair{<:Integer,<:Number}...; kwargs...) = MetricFunction(Dict(pairs); kwargs...)

"""Return the Laurent coefficients as a dictionary."""
laurent_polynomial(f::MetricFunction) = copy(f.coeffs)

"""Evaluate the metric function `f(r)`."""
function metric_value(f::MetricFunction, r::Number)
    z = ComplexF64(r)
    s = 0.0 + 0.0im
    for n in f.exponents
        s += f.coeffs[n] * z^n
    end
    return s
end

(f::MetricFunction)(r::Number) = metric_value(f, r)

"""Evaluate `df/dr` for the Laurent-polynomial metric function."""
function metric_derivative(f::MetricFunction, r::Number)
    z = ComplexF64(r)
    s = 0.0 + 0.0im
    for n in f.exponents
        n == 0 && continue
        s += n * f.coeffs[n] * z^(n - 1)
    end
    return s
end

"""
    polynomialized_metric(f)

Return `(P, m)`, where `P` is the ascending coefficient vector of

```math
P(r) = r^m f(r), \qquad m = \max(0, -n_{\min}).
```

Thus

```math
\frac{1}{f(r)} = \frac{r^m}{P(r)}.
```
"""
function polynomialized_metric(f::MetricFunction)
    m = max(0, -f.nmin)
    degree = f.nmax + m
    degree < 0 && error("polynomialized metric has negative degree")
    P = zeros(ComplexF64, degree + 1)
    for n in f.exponents
        P[n + m + 1] += f.coeffs[n]
    end
    return _trim(P), m
end

# -----------------------------------------------------------------------------
# Small polynomial toolkit, ascending coefficients.
# -----------------------------------------------------------------------------

function _trim(p::Vector{ComplexF64}; atol::Real=1.0e-14)
    q = copy(p)
    while length(q) > 1 && abs(q[end]) <= atol
        pop!(q)
    end
    return q
end

_degree(p::Vector{ComplexF64}) = length(_trim(p)) - 1

function _poly_eval(p::Vector{ComplexF64}, x::Number)
    z = ComplexF64(x)
    y = 0.0 + 0.0im
    for c in reverse(p)
        y = y * z + c
    end
    return y
end

function _poly_derivative(p::Vector{ComplexF64})
    length(p) <= 1 && return ComplexF64[0]
    q = Vector{ComplexF64}(undef, length(p) - 1)
    for k in 2:length(p)
        q[k - 1] = (k - 1) * p[k]
    end
    return _trim(q)
end

function _poly_integral_eval(p::Vector{ComplexF64}, x::Number)
    z = ComplexF64(x)
    y = 0.0 + 0.0im
    # Sum c_k z^(k+1)/(k+1), where k starts from 0.
    for k in 0:(length(p)-1)
        y += p[k + 1] * z^(k + 1) / (k + 1)
    end
    return y
end

function _poly_divrem(numer::Vector{ComplexF64}, denom::Vector{ComplexF64}; atol::Real=1.0e-12)
    n = _trim(numer; atol=atol)
    d = _trim(denom; atol=atol)
    length(d) == 1 && abs(d[1]) <= atol && error("division by zero polynomial")

    degn = length(n) - 1
    degd = length(d) - 1
    if degn < degd
        return ComplexF64[0], n
    end

    q = zeros(ComplexF64, degn - degd + 1)
    r = copy(n)
    while length(r) - 1 >= degd && !(length(r) == 1 && abs(r[1]) <= atol)
        degr = length(r) - 1
        k = degr - degd
        c = r[end] / d[end]
        q[k + 1] += c
        for j in 0:degd
            r[k + j + 1] -= c * d[j + 1]
        end
        r = _trim(r; atol=atol)
    end
    return _trim(q; atol=atol), _trim(r; atol=atol)
end

function _poly_roots(p::Vector{ComplexF64}; atol::Real=1.0e-14)
    c = _trim(p; atol=atol)
    d = length(c) - 1
    d <= 0 && return ComplexF64[]
    d == 1 && return ComplexF64[-c[1] / c[2]]

    # Normalize leading coefficient and build the companion matrix for
    # x^d + a_{d-1} x^(d-1) + ... + a_0.
    lead = c[end]
    a = c[1:end-1] ./ lead
    C = zeros(ComplexF64, d, d)
    for i in 1:(d-1)
        C[i + 1, i] = 1.0 + 0.0im
    end
    C[:, d] .= .-a
    return eigvals(C)
end

function _group_roots(roots::Vector{ComplexF64}; tol::Real=DEFAULT_ROOT_TOL)
    unused = trues(length(roots))
    groups = Vector{Vector{ComplexF64}}()
    for i in eachindex(roots)
        unused[i] || continue
        group = ComplexF64[roots[i]]
        unused[i] = false
        changed = true
        while changed
            changed = false
            center = sum(group) / length(group)
            for j in eachindex(roots)
                unused[j] || continue
                scale = max(1.0, abs(center), abs(roots[j]))
                if abs(roots[j] - center) <= tol * scale
                    push!(group, roots[j])
                    unused[j] = false
                    changed = true
                end
            end
        end
        push!(groups, group)
    end
    return groups
end

# -----------------------------------------------------------------------------
# Horizons and physical intervals.
# -----------------------------------------------------------------------------

struct Horizon
    r::Float64
    multiplicity::Int
    fp::Float64
    residue::Float64
    surface_gravity::Float64
end

struct PhysicalInterval
    left::Float64
    right::Float64
    left_kind::Symbol
    right_kind::Symbol
    sign::Int
    name::Symbol
end

Base.isfinite(I::PhysicalInterval) = isfinite(I.left) && isfinite(I.right)

struct Horizons
    roots::Vector{ComplexF64}
    groups::Vector{Vector{ComplexF64}}
    real_roots::Vector{Float64}
    positive_roots::Vector{Float64}
    horizons::Vector{Horizon}
    intervals::Vector{PhysicalInterval}
    infinity_kind::Symbol
end

"""Surface gravity at a simple horizon, `kappa = |f'(r_h)|/2`."""
surface_gravity(f::MetricFunction, r::Real) = abs(real(metric_derivative(f, r))) / 2
surface_gravity(h::Horizon) = h.surface_gravity

function _root_multiplicity(root::ComplexF64, groups::Vector{Vector{ComplexF64}}; tol::Real=DEFAULT_ROOT_TOL)
    for g in groups
        center = sum(g) / length(g)
        scale = max(1.0, abs(center), abs(root))
        if abs(root - center) <= tol * scale
            return length(g)
        end
    end
    return 1
end

function _real_roots(roots::Vector{ComplexF64}; imag_tol::Real=DEFAULT_IMAG_TOL)
    rr = Float64[]
    for z in roots
        scale = max(1.0, abs(real(z)))
        if abs(imag(z)) <= imag_tol * scale
            push!(rr, real(z))
        end
    end
    sort!(rr)
    return rr
end

function _infinity_kind(f::MetricFunction)
    # If f(r) ~ a r^p at infinity, then integral dr/f is finite at infinity
    # for p > 1.  This covers AdS-like f ~ r^2/L^2.
    return f.nmax > 1 ? :finite : :infinite
end

function _sign_real(x::ComplexF64; tol::Real=1.0e-10)
    abs(imag(x)) > tol * max(1.0, abs(real(x))) && return 0
    real(x) > 0 && return 1
    real(x) < 0 && return -1
    return 0
end

function _physical_intervals_from_roots(f::MetricFunction, positive_roots::Vector{Float64})
    intervals = PhysicalInterval[]
    roots = sort(positive_roots)

    # Intervals between neighboring positive roots.
    for i in 1:(length(roots)-1)
        a, b = roots[i], roots[i + 1]
        mid = (a + b) / 2
        s = _sign_real(metric_value(f, mid))
        name = s > 0 ? :static_finite : :nonstatic_finite
        push!(intervals, PhysicalInterval(a, b, :horizon, :horizon, s, name))
    end

    # Exterior interval beyond the largest positive root.
    if !isempty(roots)
        a = roots[end]
        probe = max(a + 1.0, 1.5a + 1.0)
        s = _sign_real(metric_value(f, probe))
        name = s > 0 ? :static_exterior : :nonstatic_exterior
        push!(intervals, PhysicalInterval(a, Inf, :horizon, :infinity, s, name))
    end

    return intervals
end

"""
    horizons(f; root_tol=1e-8, imag_tol=1e-9)

Compute roots of `P(r)=r^m f(r)`, identify positive real horizons, estimate
multiplicities, surface gravities, and candidate physical intervals.
"""
function horizons(f::MetricFunction; root_tol::Real=DEFAULT_ROOT_TOL, imag_tol::Real=DEFAULT_IMAG_TOL)
    P, _ = polynomialized_metric(f)
    roots = _poly_roots(P)
    groups = _group_roots(roots; tol=root_tol)
    rr = _real_roots(roots; imag_tol=imag_tol)
    positive = [r for r in rr if r > root_tol]

    hs = Horizon[]
    for r in positive
        mult = _root_multiplicity(ComplexF64(r), groups; tol=root_tol)
        fp = real(metric_derivative(f, r))
        residue = mult == 1 ? 1 / fp : NaN
        kappa = abs(fp) / 2
        push!(hs, Horizon(r, mult, fp, residue, kappa))
    end

    intervals = _physical_intervals_from_roots(f, positive)
    return Horizons(roots, groups, rr, positive, hs, intervals, _infinity_kind(f))
end

physical_intervals(h::Horizons; sign::Union{Nothing,Int}=nothing) =
    sign === nothing ? h.intervals : [I for I in h.intervals if I.sign == sign]

"""
    static_region(horizons; prefer=:finite)

Return a positive-`f` physical interval.  With `prefer=:finite`, this selects a
finite static interval if one exists, which is the natural Schwarzschild--dS
choice.  Otherwise it falls back to the exterior static interval.
"""
function static_region(h::Horizons; prefer::Symbol=:finite)
    candidates = physical_intervals(h; sign=1)
    isempty(candidates) && error("no positive-f static interval was found")

    finite = [I for I in candidates if isfinite(I.right)]
    exterior = [I for I in candidates if !isfinite(I.right)]

    if prefer == :finite && !isempty(finite)
        lengths = [I.right - I.left for I in finite]
        return finite[argmax(lengths)]
    elseif prefer == :exterior && !isempty(exterior)
        return exterior[end]
    else
        return candidates[end]
    end
end

static_region(f::MetricFunction; kwargs...) = static_region(horizons(f); kwargs...)

# -----------------------------------------------------------------------------
# Partial fractions and tortoise map.
# -----------------------------------------------------------------------------

struct PoleTerm
    root::ComplexF64
    order::Int
    coeff::ComplexF64
end

struct TortoiseMap
    metric::MetricFunction
    horizons::Horizons
    interval::PhysicalInterval
    terms::Vector{PoleTerm}
    quotient::Vector{ComplexF64}
    constant::ComplexF64
    convention::Symbol
end

function _partial_fraction_simple(f::MetricFunction; root_tol::Real=DEFAULT_ROOT_TOL)
    P, m = polynomialized_metric(f)
    numerator = zeros(ComplexF64, m + 1)
    numerator[m + 1] = 1.0 + 0.0im
    quotient, remainder = _poly_divrem(numerator, P)
    roots = _poly_roots(P)
    groups = _group_roots(roots; tol=root_tol)

    if any(length(g) > 1 for g in groups)
        @warn "Repeated or nearly repeated roots detected. AutoTortoise v0.1 constructs only simple-pole residues reliably. Near-extremal cases should be treated with care." max_multiplicity=maximum(length.(groups))
    end

    Pprime = _poly_derivative(P)
    terms = PoleTerm[]
    for z in roots
        denom = _poly_eval(Pprime, z)
        if abs(denom) <= 1.0e-12
            @warn "Skipping an ill-conditioned pole residue" root=z
            continue
        end
        A = _poly_eval(remainder, z) / denom
        push!(terms, PoleTerm(z, 1, A))
    end
    return terms, quotient
end

function _raw_tortoise_complex(terms::Vector{PoleTerm}, quotient::Vector{ComplexF64}, z::Number)
    zz = ComplexF64(z)
    s = _poly_integral_eval(quotient, zz)
    for t in terms
        if t.order == 1
            s += t.coeff * log(zz - t.root)
        else
            # Integral of A/(r-r0)^k for k >= 2.
            s += -t.coeff / (t.order - 1) * (zz - t.root)^(-(t.order - 1))
        end
    end
    return s
end

function _raw_tortoise_real(terms::Vector{PoleTerm}, quotient::Vector{ComplexF64}, r::Real; imag_tol::Real=DEFAULT_IMAG_TOL)
    z = ComplexF64(r)
    s = _poly_integral_eval(quotient, z)
    for t in terms
        if t.order == 1
            if abs(imag(t.root)) <= imag_tol * max(1.0, abs(real(t.root)))
                # Keep the physical real branch on a real interval.
                s += t.coeff * log(abs(float(r) - real(t.root)))
            else
                s += t.coeff * log(z - t.root)
            end
        else
            s += -t.coeff / (t.order - 1) * (z - t.root)^(-(t.order - 1))
        end
    end
    return real(s)
end

function _default_convention(f::MetricFunction, I::PhysicalInterval, h::Horizons)
    if I.right_kind == :infinity && h.infinity_kind == :finite
        return :boundary_zero
    else
        return :full_line
    end
end

"""
    tortoise_map(f, interval; convention=:auto, origin_r=nothing)

Construct the horizon-aware tortoise map

```math
r_*(r)=\int^r \frac{d\bar r}{f(\bar r)}.
```

Conventions
===========

- `:full_line`: arbitrary additive constant set to zero unless `origin_r` is
  specified.
- `:boundary_zero`: for a finite endpoint at infinity, normalize the finite
  boundary to `r_*=0`.  For standard AdS-like rational metrics this is obtained
  with zero additive constant after partial fractions.
- `:origin_zero`: set `r_*(origin_r)=0`.
"""
function tortoise_map(f::MetricFunction, I::PhysicalInterval;
                      convention::Symbol=:auto,
                      origin_r::Union{Nothing,Real}=nothing,
                      root_tol::Real=DEFAULT_ROOT_TOL)
    h = horizons(f; root_tol=root_tol)
    terms, quotient = _partial_fraction_simple(f; root_tol=root_tol)
    conv = convention == :auto ? _default_convention(f, I, h) : convention

    C = 0.0 + 0.0im
    if origin_r !== nothing || conv == :origin_zero
        origin_r === nothing && error("origin_r must be provided for convention=:origin_zero")
        C = -(_raw_tortoise_real(terms, quotient, origin_r) + 0.0im)
    elseif conv == :boundary_zero
        if !(I.right_kind == :infinity)
            error("convention=:boundary_zero requires an interval ending at infinity")
        end
        if h.infinity_kind != :finite
            error("infinity is not a finite tortoise-coordinate boundary for this metric")
        end
        # For rational functions decaying faster than 1/r, the partial-fraction
        # antiderivative tends to zero at infinity.  Keep C=0.
        C = 0.0 + 0.0im
    elseif conv == :full_line
        C = 0.0 + 0.0im
    else
        error("unknown tortoise-map convention: $conv")
    end

    return TortoiseMap(f, h, I, terms, quotient, C, conv)
end

tortoise_map(f::MetricFunction; kwargs...) = tortoise_map(f, static_region(f); kwargs...)

"""Evaluate the real-branch tortoise coordinate on the physical real interval."""
function tortoise(tm::TortoiseMap, r::Real)
    return _raw_tortoise_real(tm.terms, tm.quotient, r) + real(tm.constant)
end

"""Evaluate the complex antiderivative with the principal logarithm."""
function tortoise_complex(tm::TortoiseMap, z::Number)
    return _raw_tortoise_complex(tm.terms, tm.quotient, z) + tm.constant
end

function _left_endpoint_probe(tm::TortoiseMap, x::Real; max_shrink::Int=18)
    a = tm.interval.left
    b = tm.interval.right
    scale = isfinite(b) ? max(1.0, b - a) : max(1.0, abs(a))
    epsrel = 1.0e-10
    lo = a + epsrel * scale
    for _ in 1:max_shrink
        if tortoise(tm, lo) < x
            return lo
        end
        epsrel *= 1.0e-2
        lo = a + epsrel * scale
    end
    return lo
end

function _right_endpoint_probe_finite(tm::TortoiseMap, x::Real; max_shrink::Int=18)
    a = tm.interval.left
    b = tm.interval.right
    scale = max(1.0, b - a)
    epsrel = 1.0e-10
    hi = b - epsrel * scale
    for _ in 1:max_shrink
        if tortoise(tm, hi) > x
            return hi
        end
        epsrel *= 1.0e-2
        hi = b - epsrel * scale
    end
    return hi
end

function _right_endpoint_probe_infinite(tm::TortoiseMap, x::Real; max_grow::Int=80)
    a = tm.interval.left
    hi = max(a + 1.0, 2.0 * abs(a) + 1.0)
    for _ in 1:max_grow
        if tortoise(tm, hi) > x
            return hi
        end
        hi = a + 2.0 * (hi - a)
    end
    error("failed to bracket the inverse map on the exterior interval")
end

"""
    inverse_tortoise(tm, x; tol=1e-12, maxiter=200)

Numerically invert the real tortoise map on its physical interval using a
bracketed bisection method.  The interval is assumed to be a static region with
`f(r)>0`, so that `r_*(r)` is monotone increasing.
"""
function inverse_tortoise(tm::TortoiseMap, x::Real; tol::Real=1.0e-12, maxiter::Int=200)
    if tm.convention == :boundary_zero && tm.interval.right_kind == :infinity && x > 1.0e-12
        error("for convention=:boundary_zero the exterior coordinate usually satisfies r_* <= 0")
    end

    lo = _left_endpoint_probe(tm, x)
    hi = isfinite(tm.interval.right) ?
         _right_endpoint_probe_finite(tm, x) :
         _right_endpoint_probe_infinite(tm, x)

    xlo = tortoise(tm, lo)
    xhi = tortoise(tm, hi)
    if !(xlo <= x <= xhi)
        error("target x=$x is outside the bracketing range [$xlo, $xhi]")
    end

    for _ in 1:maxiter
        mid = (lo + hi) / 2
        xm = tortoise(tm, mid)
        if abs(xm - x) <= tol || abs(hi - lo) <= tol * max(1.0, abs(mid))
            return mid
        end
        if xm < x
            lo = mid
        else
            hi = mid
        end
    end
    return (lo + hi) / 2
end

"""
    inverse_tortoise_complex(tm, z; initial_guess=nothing)

Experimental complex inverse map.  It uses Newton iteration on the principal
logarithmic branch.  For complex-scaled CSM applications, the caller should
supply a continuous sequence of initial guesses to remain on the desired branch.
"""
function inverse_tortoise_complex(tm::TortoiseMap, z::Number;
                                  initial_guess::Union{Nothing,Number}=nothing,
                                  tol::Real=1.0e-12,
                                  maxiter::Int=60)
    target = ComplexF64(z)
    r0 = if initial_guess === nothing
        try
            inverse_tortoise(tm, real(target)) + 0.0im
        catch
            if isfinite(tm.interval.right)
                (tm.interval.left + tm.interval.right) / 2 + 0.0im
            else
                tm.interval.left + 1.0 + 0.0im
            end
        end
    else
        ComplexF64(initial_guess)
    end

    r = r0
    for _ in 1:maxiter
        F = tortoise_complex(tm, r) - target
        abs(F) <= tol && return r
        # d r_*/d r = 1/f(r), hence Newton update is r <- r - F*f(r).
        r -= F * metric_value(tm.metric, r)
    end
    error("complex inverse_tortoise did not converge; try a better initial_guess or branch continuation")
end

"""
    potential_at(tm, V, x; theta=0, initial_guess=nothing)

Evaluate a potential `V(r)` at `r = r(r_* exp(i theta))`.  This is a small
CSM-oriented convenience wrapper around `inverse_tortoise_complex`.
"""
function potential_at(tm::TortoiseMap, V::Function, x::Real;
                      theta::Real=0.0,
                      initial_guess::Union{Nothing,Number}=nothing)
    zstar = ComplexF64(x) * exp(1im * theta)
    r = inverse_tortoise_complex(tm, zstar; initial_guess=initial_guess)
    return V(r)
end

# -----------------------------------------------------------------------------
# Common metric constructors.
# -----------------------------------------------------------------------------

"""Schwarzschild metric function, `f(r)=1-2M/r`."""
schwarzschild(; M::Real=1.0) = MetricFunction(-1 => -2M, 0 => 1.0)

"""Schwarzschild with cosmological constant, `f(r)=1-2M/r-Λ r^2/3`."""
schwarzschild_lambda(; M::Real=1.0, Lambda::Real=0.0) =
    MetricFunction(-1 => -2M, 0 => 1.0, 2 => -Lambda / 3)

"""Schwarzschild--de Sitter metric function with positive `Lambda`."""
schwarzschild_ds(; M::Real=1.0, Lambda::Real=0.01) =
    schwarzschild_lambda(; M=M, Lambda=Lambda)

"""
Schwarzschild--AdS metric function.  Either pass `L`, giving
`f(r)=1-2M/r+r^2/L^2`, or pass a negative `Lambda`.
"""
function schwarzschild_ads(; M::Real=1.0, L::Union{Nothing,Real}=nothing, Lambda::Union{Nothing,Real}=nothing)
    if L !== nothing
        return MetricFunction(-1 => -2M, 0 => 1.0, 2 => 1 / L^2)
    elseif Lambda !== nothing
        Lambda >= 0 && error("Schwarzschild--AdS requires negative Lambda")
        return schwarzschild_lambda(; M=M, Lambda=Lambda)
    else
        error("provide either L or negative Lambda")
    end
end

end # module
