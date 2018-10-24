module GenomicMaps

using Luxor
using Colors
using GenomicAnnotations
using Statistics
using Printf: @sprintf

export drawgenome

include("shapes.jl")
include("drawstretch.jl")

end #module
