module PhaseSectors

export PhaseSector, 𝕋, dist, unwrap, dilate, sectors!

include("PhaseSector.jl")

include("helpers.jl")
include("utils.jl")

function sectors! #= ext/PhaseSectorsMakieExt.jl =# end

end # module PhaseSectors
