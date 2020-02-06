@doc raw"""
Contains functions and types to work with the analytic solutions to the Dirac equation.

Defines the functions `f` and `g` that calculate the radial functions of the solution of
the single-particle Dirac equation in a ``1/r`` potential. The `f` and `g` are defined as

```math
\psi(r,\theta,\varphi)
= \frac{1}{r}
\begin{pmatrix}
g(r) \chi_{\kappa m}(\theta,\varphi) \\
i f(r) \chi_{-\kappa m}(\theta, \varphi)
\end{pmatrix}
```

Atomic units, i.e. ``c = 1/α``, ``μ = 1`` and ``\hbar = 1``.

<https://en.wikipedia.org/wiki/Hydrogen-like_atom#Solution_to_Dirac_equation>
"""
module DiracAnalytic
using AtomicLevels: RelativisticOrbital

"""
The fine structure constants ``\\alpha``.
"""
const α = 0.0072973525664

@doc raw"""
    E(Z, n, κ) -> Real
    E(Z, orbital::RelativisticOrbital) -> Real

Calculates the relativistic energy of a hydrogenic orbital in the ``V(r) = - Z / \alpha``
potential:

```math
E_{n\kappa}(Z) =
m c^2\left(1+\left[\frac{Z\alpha}{n-|k|+\sqrt{k^2-Z^2\alpha^2}}\right]^2\right)^{-1/2}
```

Also supports the following signatures

    E(::Type{Complex}, ...) -> Complex

which allow the value to become complex. This happens after the critical nuclear charge
``Z > \abs{\kappa} / \alpha``.

```julia-repl
julia> E(Complex, 200, ro"1s")
0.0 + 19962.685927945768im
```
"""
function E end
E(Z::Real, orb::RelativisticOrbital) = E(Z, orb.n, orb.κ)
E(T::Type{Complex}, Z, orb::RelativisticOrbital) = E(T, Z, orb.n, orb.κ)
E(Z::Real, n::Integer, κ::Integer) = _E(Z, n, κ, _Δ(Z, κ))
E(::Type{Complex}, Z, n::Integer, κ::Integer) = _E(Z, n, κ, complex(_Δ(Z, κ)))
# Internal methods to handle the complex and non-complex versions of E
_E(Z, n, κ, Δ) = 1 / (α^2 * sqrt(1 + (Z*α / (n - abs(κ) + sqrt(Δ)))^2))
_Δ(Z, κ) = κ^2 - (Z*α)^2

end # module
