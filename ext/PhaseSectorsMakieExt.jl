module PhaseSectorsMakieExt

using PhaseSectors

using Makie: Axis, Colorant, lines!, band!, translate!

include("helpers.jl")

"""
    sectors!(ax, x, Φ[; linecolor, fillcolor])

Plot vector `Φ` of phase sectors against vector `x` on Makie axis `ax`. The sectors are
unwrapped internally, and the plot is repeated with 2π shifts, so the vertical axis can
be cropped to [-π, π].
"""
function PhaseSectors.sectors!(
        ax::Axis, x::AbstractVector{<:Real}, Φ::AbstractVector{<:PhaseSector};
        linecolor = :black, fillcolor = (:black, 0.2)
    )
    Φun = unwrap(Φ)
    αun, βun = [getproperty.(Φun, i) for i in (:α, :β)]
    αdr, βdr = [getproperty.(_drawable(Φun), i) for i in 1:2]
    lowest_band = minimum(floor.((αdr .+ π) / 2π))
    highest_band = maximum(ceil.((βdr .- π) / 2π))
    for shift in 2π * ((-highest_band - 1):(-lowest_band + 1))
        p = band!(ax, x, αdr, βdr; color = fillcolor)
        p = [p; lines!(ax, x, αun; color = linecolor)]
        p = [p; lines!(ax, x, βun; color = linecolor)]
        @views for (x, α, β, Φprev, Φcur, Φnext) in zip(
                x[2:end], αdr[2:end], βdr[2:end],
                Φun[1:end], Φun[2:end], Φun[3:end]
            )
            if Φcur ≠ 𝕋 && (Φprev ≡ 𝕋 || Φnext ≡ 𝕋)
                p = [p; lines!(ax, [x, x], [β, α + 2π]; color = linecolor)]
            end
        end
        translate!.(p, 0, shift, 0)
    end
    return nothing
end

end # module PhaseSectorsMakieExt
