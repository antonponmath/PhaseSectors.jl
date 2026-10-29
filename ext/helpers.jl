_drawable(Φ::AbstractVector{PhaseSector}) = accumulate(
    (prev, cur) ->
    cur ≡ 𝕋 ? (prev[1], prev[1] + 2π) : (cur.α, cur.β),
    Φ; init = (0, 2π)
)
