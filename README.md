# PhaseSectors.jl

Operations with phase sectors as defined in [[1]](#1). Intended for the small-phase analysis of control systems.

## Definitions

**Phases** are elements of the space $\mathbb{T} = \mathbb{R}/(2\pi\mathbb{Z})$. Phase $\theta \in \mathbb{T}$ is therefore an equivalence class. Notation $\theta = [x]$ means that $x \in \mathbb{R}$ is a representative of that class.

**Phase arc** covered by an interval $[x_1, x_2] \subset \mathbb{R}$ is
```math
\mathrm{arc}[x_1, x_2] = \bigl\{ [x] \colon x \in [x_1, x_2] \bigr\}.
```

**Metric** on $\mathbb{T}$ is
```math
\mathrm{dist}(\theta_1, \theta_2) = \min\bigl\{ |x_1 - x_2| \colon [x_1] = \theta_1, [x_2] = \theta_2 \bigr\}.
```

The **diameter** of a set $\Theta \subseteq \mathbb{T}$ is
```math
\mathrm{diam}\Theta = \sup\bigl\{ \mathrm{dist}(\theta_1, \theta_2) \colon \theta_1, \theta_2 \in \Theta \bigr\}.
```

Set $\Theta \subseteq \mathbb{T}$ is called **convex** if every pair of its points can be connected by a minimizing $\mathbb{T}$-geodesic contained in $\Theta$.

**Phase sector** $\Phi$ is a convex phase arc: that is,
- either the entire circle: $\Phi = \mathbb{T}$,
- or a phase arc with $\mathrm{diam}\Phi < \pi$.

## What Is Implemented

Create the phase sector $\Phi = \mathrm{arc}[\alpha, \beta]$:
```julia
Φ = PhaseSector(α, β)
```
If $\beta \geq \alpha + \pi$, then the sector automatically becomes the entire circle $\mathbb{T}$ which is defined as
```julia
const 𝕋 = PhaseSector(-Inf, Inf)
```

Endpoints `Φ.α` and `Φ.β` can be also accessed as `Φ.a` and `Φ.b`.

**Operations:**
- Check if a phase is inside the sector: `θ ∈ Φ`.
- Minkowski sum of phase sectors, etc.: `Φ1 + Φ2`, `Φ + θ`, `Φ1 - Φ2`, ...
- Diameter of a phase sector: `diam(Φ)`.
- Distance between a phase and a sector, or between sectors: `dist(θ, Φ)`, `dist(Φ1, Φ2)`.
- Unwrap a sequence of phase sectors similarly to `unwrap` for angles but ignoring entire $\mathbb{T}$-sectors: `unwrap(Φ)`. Useful for plotting.
- Join (convex hull of the union) of sectors: `join([Φ1, Φ2, ...])` or `join(Φ1, Φ2, ...)`.
- Join each sector in a sequence with its neighbors: `dilate([Φ1, Φ2, ...])`. This makes sure that neighboring sectors overlap, and thus the sequence is "continuous." In small-phase analysis, it is helpful for zero crossing detection: find out if $0 \in \Phi_k$ for some $k$.
- If `Makie` is loaded, then the plotting method `sector!` is available to plot sequences of phase sectors.

---

<a id="1">[1]</a> A. Ponomarev. *Compositional small-gain and small-phase stability analysis*. 2026. [arXiv:2606.14745](https://arxiv.org/abs/2606.14745).
