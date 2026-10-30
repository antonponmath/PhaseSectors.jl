"""
    PhaseSector

Phase sector described by its endpoint angles `α` and `β` or, equivalently, `a` and `b`.

    PhaseSector{T}(α::Real, β::Real) where {T <: Real}
    PhaseSector(α::Real, β::Real)
    PhaseSector(αβ::Tuple{<:Real, <:Real})

Construct the phase sector from its endpoint angles `α` and `β ≥ α`, or their tuple
`αβ`. If `β ≥ α + π`, the sector becomes infinite (the entire 2π circle, denoted `𝕋`).

"""
struct PhaseSector{T <: Real}
    α::T
    β::T
    function PhaseSector{T}(α::Real, β::Real) where {T <: Real}
        α, β = T(α), T(β)
        α ≤ β || throw(ArgumentError("phase bounds must be ordered, got ($α, $β)"))
        if β ≥ α + π
            α, β = -Inf, Inf
        end
        return new(α, β)
    end
end

PhaseSector(α::T, β::T) where {T <: Real} = PhaseSector{T}(α, β)

PhaseSector(α::Real, β::Real) = PhaseSector(promote(α, β)...)

PhaseSector(αβ::Tuple{<:Real, <:Real}) = PhaseSector(αβ...)


"""
    PhaseSector(z::AbstractVector{<:Number})

Construct the phase sector spanned by a set of complex numbers.
"""
PhaseSector(z::AbstractVector{<:Number}) =
    z .|> angle |> _unwrap |> extrema |> PhaseSector


"""
    𝕋 ≡ PhaseSector(-Inf, Inf)

The infinite phase sector representing the entire 2π circle.
"""
const 𝕋 = PhaseSector(-Inf, Inf)


Base.broadcastable(Φ::PhaseSector) = Ref(Φ)


Base.print(io::IO, Φ::PhaseSector) = begin
    if Φ ≡ 𝕋
        print(io, "𝕋")
    else
        Φ = Φ - _div2pi((Φ.α + Φ.β) / 2)
        print(io, "arc[$(Φ.α / π)π, $(Φ.β / π)π]")
    end
end


"Defining aliases a and b for α and β."
Base.getproperty(Φ::PhaseSector, name::Symbol) =
    name ≡ :a ? getfield(Φ, :α) :
    name ≡ :b ? getfield(Φ, :β) : getfield(Φ, name)


Base.:-(Φ::PhaseSector) = PhaseSector(-Φ.β, -Φ.α)

Base.:+(Φ::PhaseSector, θ::Real) = PhaseSector(Φ.α + θ, Φ.β + θ)
Base.:-(Φ::PhaseSector, θ::Real) = Φ + (-θ)
Base.:+(θ::Real, Φ::PhaseSector) = Φ + θ
Base.:-(θ::Real, Φ::PhaseSector) = (-Φ) + θ

Base.:+(Φ1::PhaseSector, Φ2::PhaseSector) = PhaseSector(Φ1.α + Φ2.α, Φ1.β + Φ2.β)
Base.:-(Φ1::PhaseSector, Φ2::PhaseSector) = Φ1 + (-Φ2)


Base.in(θ::Real, Φ::PhaseSector) =
    Φ ≡ 𝕋 || (mod2pi(θ - Φ.α) < π && mod2pi(Φ.β - θ) < π)
