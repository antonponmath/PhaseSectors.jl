module PhaseSectors

export PhaseSector, 𝕋, dist, boundary, unwrap, dilate


include("helpers.jl")

struct PhaseSector
    α::Real
    β::Real
    a::Real # ≡ α
    b::Real # ≡ β
    function PhaseSector(α::Real, β::Real)
        α ≤ β || throw(ArgumentError("phase bounds must be ordered, got ($α, $β)"))
        if β ≥ α + π
            α, β = -Inf, Inf
        end
        return new(α, β, α, β)
    end
end


𝕋 = PhaseSector(-Inf, Inf)


PhaseSector(αβ::Tuple{<:Real, <:Real}) = PhaseSector(αβ...)


PhaseSector(z::AbstractVector{<:Number}) =
    z .|> angle |> _unwrap |> extrema |> PhaseSector


Base.broadcastable(Φ::PhaseSector) = Ref(Φ)


Base.show(io::IO, Φ::PhaseSector) =
    print(io, "PhaseSector($(Φ.α), $(Φ.β))")


Base.print(io::IO, Φ::PhaseSector) = begin
    Φ = Φ - _div2pi((Φ.α + Φ.β) / 2)
    print(io, Φ ≡ 𝕋 ? "𝕋" : "arc[$(Φ.α / π)π, $(Φ.β / π)π]")
end


Base.:-(Φ::PhaseSector) = PhaseSector(-Φ.β, -Φ.α)

Base.:+(Φ::PhaseSector, θ::Real) = PhaseSector(Φ.α + θ, Φ.β + θ)
Base.:-(Φ::PhaseSector, θ::Real) = Φ + (-θ)
Base.:+(θ::Real, Φ::PhaseSector) = Φ + θ
Base.:-(θ::Real, Φ::PhaseSector) = (-Φ) + θ

Base.:+(Φ1::PhaseSector, Φ2::PhaseSector) = PhaseSector(Φ1.α + Φ2.α, Φ1.β + Φ2.β)
Base.:-(Φ1::PhaseSector, Φ2::PhaseSector) = Φ1 + (-Φ2)


Base.in(θ::Real, Φ::PhaseSector) =
    Φ ≡ 𝕋 || (mod2pi(θ - Φ.α) < π && mod2pi(Φ.β - θ) < π)


function dist(θ::Number, Φ::PhaseSector)
    θ ∈ Φ && return 0
    return min(mod2pi(Φ.α - θ), mod2pi(θ - Φ.β))
end


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


boundary(Φ::PhaseSector) = (Φ.α, Φ.β)


function unwrap(Φ::AbstractVector{PhaseSector})
    Φun = copy(Φ)
    mask = (Φ .≠ 𝕋)
    Φun[mask] .= accumulate((prev, cur) -> cur - _div2pi(cur.α - prev.α), Φ[mask])
    return Φun
end


function Base.join(Φ::AbstractVector{PhaseSector})
    Φun = unwrap(Φ)
    return PhaseSector(minimum(getproperty.(Φun, :α)), maximum(getproperty.(Φun, :β)))
end


Base.join(Φ1::PhaseSector, Φ2::PhaseSector) = join([Φ1, Φ2])


dilate(Φ::AbstractVector{PhaseSector}) = [
    join(@view Φ[((i - 1):(i + 1)) ∩ eachindex(Φ)]) for i in eachindex(Φ)
]

end # module PhaseSectors
