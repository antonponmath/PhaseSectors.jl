module PhaseSectorsMakieExt

using PhaseSectors

using Makie: Axis, lines!, band!, translate!


_drawable(Φ::AbstractVector{PhaseSector}) = accumulate(
    (prev, cur) ->
    cur ≡ 𝕋 ? (prev[1], prev[1] + 2π) : (cur.α, cur.β),
    Φ; init = (0, 2π)
)

function PhaseSectors.sectors!(
        ax::Axis, ω::AbstractVector{<:Real}, Φ::AbstractVector{PhaseSector};
        linecolor = :black, fillcolor = (:black, 0.2)
    )
    Φun = unwrap(Φ)
    αun, βun = [getproperty.(Φun, i) for i in (:α, :β)]
    αdr, βdr = [getproperty.(_drawable(Φun), i) for i in 1:2]
    lowest_band = minimum(floor.((αdr .+ π) / 2π))
    highest_band = maximum(ceil.((βdr .- π) / 2π))
    for shift in 2π * ((-highest_band - 1):(-lowest_band + 1))
        p = band!(ax, ω, αdr, βdr; color = fillcolor)
        p = [p; lines!(ax, ω, αun; color = linecolor)]
        p = [p; lines!(ax, ω, βun; color = linecolor)]
        @views for (ω, α, β, Φprev, Φcur, Φnext) in zip(
                ω[2:end], αdr[2:end], βdr[2:end],
                Φun[1:end], Φun[2:end], Φun[3:end]
            )
            if Φcur ≠ 𝕋 && (Φprev ≡ 𝕋 || Φnext ≡ 𝕋)
                p = [p; lines!(ax, [ω, ω], [β, α + 2π]; color = linecolor)]
            end
        end
        translate!.(p, 0, shift, 0)
    end
    return
end

end # module PhaseSectorsMakieExt
