"""
    dist(θ::Real, Φ::PhaseSector)
    dist(Φ::PhaseSector, θ::Real)

Distance from an angle to a phase sector in radians.
"""
function dist(θ::Real, Φ::PhaseSector)
    θ ∈ Φ && return 0
    return min(mod2pi(Φ.α - θ), mod2pi(θ - Φ.β))
end

dist(Φ::PhaseSector, θ::Real) = dist(θ, Φ)

"""
    dist(Φ1::PhaseSector, Φ2::PhaseSector)

Distance between two phase sectors in radians.
"""
function dist(Φ1::PhaseSector, Φ2::PhaseSector)
    (Φ1 ≡ 𝕋 || Φ2 ≡ 𝕋) && return 0
    return minimum(
        map(
            dist,
            [Φ1.α, Φ1.β, Φ2.α, Φ2.β],
            [Φ2, Φ2, Φ1, Φ1]
        )
    )
end


"""
    unwrap(Φ)

Unwrap a sequence of phase sectors similarly to unwrapping angles. Ignores infinite
sectors.
"""
function unwrap(Φ::AbstractVector{PhaseSector})
    Φun = copy(Φ)
    mask = (Φ .≠ 𝕋)
    Φun[mask] .= accumulate(
        (prev, cur) -> cur - _div2pi(cur.α - prev.α),
        Φ[mask]
    )
    return Φun
end


"""
    join(Φ::AbstractVector{PhaseSector})
    join(Φ1::PhaseSector, Φ2::PhaseSector, etc::PhaseSector...)

The join of phase sectors (convex hull of their union).
"""
function Base.join(Φ::AbstractVector{PhaseSector})
    Φun = unwrap(Φ)
    return PhaseSector(
        minimum(getproperty.(Φun, :α)),
        maximum(getproperty.(Φun, :β))
    )
end

Base.join(Φ1::PhaseSector, Φ2::PhaseSector, etc::PhaseSector...) =
    join(PhaseSector[Φ1, Φ2, etc...])


"""
    dilate(Φ::AbstractVector{PhaseSector})

Join each phase sector with its neighbors. Helps detecting zero crossings when the
sectors are thin and the zero can slip between them.
"""
dilate(Φ::AbstractVector{PhaseSector}) = [
    join(@view Φ[((i - 1):(i + 1)) ∩ eachindex(Φ)]) for i in eachindex(Φ)
]
