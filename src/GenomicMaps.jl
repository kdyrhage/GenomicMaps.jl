module GenomicMaps

using Luxor
using Colors
using GenomicAnnotations
using GenomicAnnotations: AbstractGene, Gene
using Statistics
using Printf: @sprintf


include("initialise.jl")
include("shapes.jl")
include("drawstretch.jl")
export drawgenome
include("annotatesvg.jl")
include("drawplasmid.jl")
export drawplasmid
include("genomecomparison.jl")
export drawcomparison

end #module
