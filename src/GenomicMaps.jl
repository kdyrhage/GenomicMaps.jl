module GenomicMaps

using Luxor
using Colors
using GenomicAnnotations
using GenomicAnnotations: AbstractGene, Gene
using Statistics
using Printf: @sprintf

export drawgenome

include("initialise.jl")
include("shapes.jl")
include("drawstretch.jl")
include("annotatesvg.jl")

end #module
