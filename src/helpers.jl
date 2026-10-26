_div2pi(θ::Real) = θ - rem2pi(θ, RoundNearest)

_unwrap(θ::AbstractVector{<:Real}) =
    accumulate((prev, cur) -> cur - _div2pi(cur - prev), θ)
