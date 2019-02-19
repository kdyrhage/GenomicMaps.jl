# GenomicMaps.jl
Julia package for visualising genomic data.

## Example
Open the image in its own tab, then hover over a gene to see annotations.
# <img src="./example.svg" />

```julia
using GenomicAnnotations
using GenomicMaps


# You can add any kind of annotation that you want to display.
# Here I add the results from InterProScan:
function addinterproscan!(chr, filename)
    for line in split.(readlines(filename), Ref('\t'))
        g = first(chr.genes[(chr.genedata[:feature] .== "CDS") .& (chr.genedata[:locus_tag] .== line[1])])
        if length(line) < 13 && isempty(line[6])
            pushproperty!(g, Symbol(line[4]), String(line[5]))
        elseif length(line) < 13
            pushproperty!(g, Symbol(line[4]), "$(line[5]): $(line[6])")
        else
            pushproperty!(g, Symbol(line[4]), "$(line[5]): $(line[13])")
        end
    end
    delete!(chr.genedata, :SUPERFAMILY)
    delete!(chr.genedata, :MobiDBLite)
    delete!(chr.genedata, :Gene3D)
end


# ... and COGs:
function addcogs!(chr, filename)
    chr.genes.cog .= ""
    @progress for line in split.(readlines(filename), Ref('\t'))
        g = first(@genes(chr, :feature == "CDS", :locus_tag == line[1]))
        if occursin(r"\w", line[2])
            g.cog = replace(line[2], r"\W" => s"")
        end
    end
    @genes(chr, isempty(:cog)).cog .= missing
end


chr = readgbk("GCA_000005845.2_ASM584v2_genomic.gbff.gz")[1]
addinterproscan!(chr, "interproresults.tsv")
addcogs!(chr, "cogresults.tsv")


# Colour scheme for COG categories:
cogcolours = Dict("B"=>RGB{Float64}(1.0,0.630714,0.576563),
     "M"=>RGB{Float64}(0.756869,0.916499,0.965176),
     "I"=>RGB{Float64}(0.187839,0.54561,0.252343),
     "X"=>RGB{Float64}(0.540006,0.493982,0.813567),
     "Y"=>RGB{Float64}(0.0973617,0.285282,0.5329),
     "Z"=>RGB{Float64}(0.0418427,0.156645,0.341597),
     "L"=>RGB{Float64}(0.426131,0.0441442,0.0465628),
     "O"=>RGB{Float64}(0.518954,0.802339,0.930272),
     "F"=>RGB{Float64}(0.587882,0.865532,0.51112),
     "Q"=>RGB{Float64}(0.0,0.225356,0.101282),
     "D"=>RGB{Float64}(0.862653,0.958477,0.981395),
     "V"=>RGB{Float64}(0.188382,0.529206,0.795898),
     "U"=>RGB{Float64}(0.277786,0.635283,0.863472),
     "E"=>RGB{Float64}(0.711814,0.932724,0.629136),
     "T"=>RGB{Float64}(0.394211,0.72627,0.90426),
     "H"=>RGB{Float64}(0.32729,0.673206,0.326717),
     "P"=>RGB{Float64}(0.0232916,0.395886,0.180144),
     "G"=>RGB{Float64}(0.459895,0.779462,0.41097),
     "N"=>RGB{Float64}(0.641543,0.865092,0.94902),
     "K"=>RGB{Float64}(0.825431,0.118066,0.106858),
     "C"=>RGB{Float64}(0.835916,0.980813,0.770886),
     "R"=>RGB{Float64}(0.807625,0.787968,0.949453),
     "W"=>RGB{Float64}(0.137797,0.411028,0.686187),
     "A"=>RGB{Float64}(1.0,0.808314,0.771835),
     "S"=>RGB{Float64}(0.236943,0.0166779,0.407047),
     "J"=>RGB{Float64}(1.0,0.389569,0.336934))

# The output can be customised, see src/initialise.jl for all options. Here I
# provide a function that will be run on each gene to determine its colour:
colourby_cog = g -> unique(String.(split(get(g, :cog, ""), "")))
drawgenome(chr;
    outfile = "example.svg",
    colourmap = cogcolours,
    colourfunction = colourby_cog,
    annotate = true,
    nbreaks = 40)
```

## Customisation
Some keywords that can be given to `drawgenome` to customise the output are:
- `genetextfunction`: determines the text shown above each gene. Can be either a `Function` that is executed for each gene, or a `Symbol`, in which case it defaults to `g -> get(g, genetextfunction, "")`.
- `colourfunction`: a function that is executed for each gene to determine how to colour it. Currently only categorical colouring is supported, so continuous data has to be binned.
- `colourmap`: a `Dict` mapping categories => colours. When using categorical data to colour genes, it can be left empty, but if continuous data (say, expression levels) are used it has to be set manually.
- `defaultcolour`: the default colour that is used for genes that do not have a specified colour in `colourmap`.
- `nbreaks`: an `Int` determining how many lines will be used to display the genome).
- `drawingsize`: determines the size of the output. Can be a `String` such as `"A4"`, `"A0landscape"`, `"1000x1000"`, or a `Tuple` (e.g. `(1000, 1000)`).
... and more (see `src/initialise.jl`)
