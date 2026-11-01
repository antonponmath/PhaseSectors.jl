"""
    diam(Φ)

Return the diameter of phase sector Φ.
"""
diam(Φ::PhaseSector) =
    return Φ ≡ 𝕋 ? π : Φ.β - Φ.α


"""
    dist(θ, Φ)
    dist(Φ, θ)

Return the distance from angle `θ` to phase sector `Φ`.
"""
function dist(θ::Real, Φ::PhaseSector)
    θ ∈ Φ && return 0
    return min(mod2pi(Φ.α - θ), mod2pi(θ - Φ.β))
end

dist(Φ::PhaseSector, θ::Real) = dist(θ, Φ)

"""
    dist(Φ1, Φ2)

Return the distance between phase sectors `Φ1` and `Φ2`.
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

Unwrap and return the vector `Φ` of phase sectors. Works similarly to unwrapping angles.
Ignores infinite sectors.
"""
function unwrap(Φ::AbstractVector{<:PhaseSector})
    Φun = copy(Φ)
    mask = (Φ .≠ 𝕋)
    Φun[mask] .= accumulate(
        (prev, cur) -> cur - _div2pi(cur.α - prev.α),
        Φ[mask]
    )
    return Φun
end


"""
    join(Φ)

Return the join of phase sectors contained in vector `Φ`.
"""
function Base.join(Φ::AbstractVector{<:PhaseSector})
    Φun = unwrap(Φ)
    return PhaseSector(
        minimum(getproperty.(Φun, :α)),
        maximum(getproperty.(Φun, :β))
    )
end

"""
    join(Φ1, Φ2[, ...])

Return the join of phase sectors `Φ1`, `Φ2`, ...
"""
Base.join(Φ1::PhaseSector, Φ2::PhaseSector, etc::PhaseSector...) =
    join(PhaseSector[Φ1, Φ2, etc...])


"""
    dilate(Φ)

Join each phase sector in vector `Φ` with its immediate neighbors and return the result.
"""
dilate(Φ::AbstractVector{<:PhaseSector}) = [
    join(@view Φ[((i - 1):(i + 1)) ∩ eachindex(Φ)]) for i in eachindex(Φ)
]
