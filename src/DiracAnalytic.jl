@doc raw"""
Contains functions and types to work with the analytic solutions to the Dirac equation.

The bound state solutions to the Dirac equation are represented via the ``P_{n\kappa}(r)``
and ``Q_{n\kappa}(r)`` functions defined for each ``n`` and ``\kappa`` quantum number
combination and defined via the relation:

```math
\psi_{n\kappa m}(r,\theta,\varphi)
= \frac{1}{r} \begin{pmatrix}
  P_{n\kappa}(r) \chi_{\kappa m}(\theta,\varphi) \\
  i Q_{n\kappa}(r) \chi_{-\kappa m}(\theta, \varphi)
\end{pmatrix}
```

!!! note "Atomic units"

    All values are in **atomic units**, i.e. ``c = 1/α``, ``m = 1`` and ``\hbar = 1``.
"""
module DiracAnalytic
using SpecialFunctions: gamma
using GSL: sf_laguerre_n
using AtomicLevels: RelativisticOrbital

"""
The fine structure constants ``\\alpha``.

```@repl
α
1/α
```
"""
const α = 0.0072973525664

@doc raw"""
    energy(Z, n, κ) -> Real
    energy(Z, orbital::RelativisticOrbital) -> Real

Calculates the relativistic energy of a hydrogenic orbital in the ``V(r) = - Z / \alpha``
potential:

```math
E_{n\kappa}(Z) =
m c^2\left(1+\left[\frac{Z\alpha}{n-|k|+\sqrt{k^2-Z^2\alpha^2}}\right]^2\right)^{-1/2}
```

Also supports the following signatures

    energy(::Type{Complex}, ...) -> Complex

which allow the value to become complex. This happens after the critical nuclear charge
``Z > \abs{\kappa} / \alpha``.

```julia-repl
julia> energy(Complex, 200, ro"1s")
0.0 + 19962.685927945768im
```
"""
function energy end
energy(Z::Real, orb::RelativisticOrbital) = energy(Z, orb.n, orb.κ)
energy(T::Type{Complex}, Z, orb::RelativisticOrbital) = energy(T, Z, orb.n, orb.κ)
energy(Z::Real, n::Integer, κ::Integer) = _E(Z, n, κ, _Δ(Z, κ))
energy(::Type{Complex}, Z, n::Integer, κ::Integer) = _E(Z, n, κ, complex(_Δ(Z, κ)))
# Internal methods to handle the complex and non-complex versions of E
_E(Z, n, κ, Δ) = 1 / (α^2 * sqrt(1 + (Z*α / (n - abs(κ) + sqrt(Δ)))^2))
_Δ(Z, κ) = κ^2 - (Z*α)^2

# Radial solutions of the DiracEquation

abstract type DiracSolution end
function PQ end
Base.broadcastable(x::DiracSolution) = Ref(x)

"""
Implementation based on the expressions in Wikipedia:
<https://en.wikipedia.org/wiki/Hydrogen-like_atom#Solution_to_Dirac_equation>
"""
struct WikiSolution <: DiracSolution
    Z :: Float64
    n :: Int
    κ :: Int
end

function PQ(s::WikiSolution, r)
    E = energy(s.Z, s.n, s.κ)
    C, γ = _wiki_C(s.Z, s.n), _wiki_γ(s.Z, s.κ)
    A = _wiki_A(s.Z, s.n, s.κ)
    ρ = 2 * C * r

    if s.n == -s.κ
        A * (s.n + γ) * (ρ^γ) * exp(-ρ/2), A * s.Z * α * (ρ^γ) * exp(-ρ/2)
    else
        L1 = sf_laguerre_n(s.n - abs(s.κ) - 1, 2γ + 1, ρ)
        L2 = sf_laguerre_n(s.n - abs(s.κ), 2γ + 1, ρ)
        A * (ρ^γ) * exp(-ρ/2) * (
            s.Z * α * ρ * L1 +
            (γ - s.κ) * (γ/α^2 - s.κ*E) * L2 * α / C
        ),
        A * (ρ^γ) * exp(-ρ/2) * (
            (γ - s.κ) * ρ * L1 +
            s.Z * α * (γ/α^2 - s.κ*E) * L2 * α / C
        )
    end
end

function _wiki_A(Z, n, k)
    C, γ = _wiki_C(Z, n), _wiki_γ(Z, k)
    if n == -k
        sqrt(C / (2 * n * (n + γ) * γ * gamma(2γ)))
    else
        E = energy(Z, n, k)
        sqrt(
            0.5 *
            (C / (n - abs(k) + γ)) *
            (factorial(n - abs(k) - 1) / gamma(n - abs(k) + 2γ + 1)) *
            ((α^2 * E * k / γ)^2 + (α^2 * E * k / γ))
        ) / sqrt(2k * (k - γ))
    end
end

@doc raw"""
    _wiki_γ(Z, κ)

Calculates the ``\gamma = \sqrt{\kappa^2 - Z^2 \alpha^2}`` value.
"""
_wiki_γ(Z, κ) = sqrt(κ^2 - (Z*α)^2)

@doc raw"""
    _wiki_C(Z, n)

Calculates the radial decay constant ``C``:

```math
C = \frac{Z\alpha}{n}\frac{m c^2}{\hbar c} = \frac{Z}{n}
```
"""
_wiki_C(Z, n) = Z/n # μ = 1, ħ = 1 and α*c^2 / c = 1

end # module
