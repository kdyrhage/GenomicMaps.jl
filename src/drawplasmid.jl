function geneangles(gene::AbstractGene, chrlength::Int, ORI::Int = 1)
    loc = locus(gene)
    geneangles(loc.position.start, loc.position.stop, chrlength, ORI)
end
function geneangles(start::Number, stop::Number, chrlength::Int, ORI::Int = 1)
    a1 = mod2pi((2π * (start-ORI-1) / chrlength)) - π/2
    a2 = mod2pi((2π * (stop-ORI-1) / chrlength)) - π/2
    amid = mod2pi((2π * (mean([start, stop])-ORI-1) / chrlength)) - π/2
    a1, a2, amid
end

"""
    function drawplasmid(chr; kwargs...)

Draw a circular plasmid map. CDSs are displayed as boxes on the outside and inside of the plasmid, for genes on the + and - strands, respectively.
Supported kwargs are:
* outfile           Path to output file. Supported file types are PNG, SVG, PDF, and EPS.
* drawingsize       Size of the figure, can be a `Tuple{Int,Int}`, or a `String` such as "1000x1000".
* title             Name of the plasmid.
* ORI               The position of the ORI (defaults to 1).
* colors            A `Dict`, with `Gene`s as keys and colors as values (defaults to `defaultgenecolor`).
* defaultgenecolor  Defaults to :lightgrey.
* annotations       A `Vector` containing annotations for each CDS.
* internalcolors    A `Dict` with `UnitRange`s (representing nucleotide positions) as keys, and colors as values. Can be used to color parts of genes, or intergenic regions.
* genethickness     Size of the gene boxes.
* textoffset        Distance from the plasmid where annotations should be drawn.
* figureoffset      A `Tuple{Int,Int}` with offsets in the x and y axes, respectively. Useful if there are more annotations on one side than the other.
* genegroups        A `Dict` with `UnitRanges` as keys (representing gene numbers), and `NamedTuples` as values. The tuple must have the field `text`, which contains the annotations for the group, and can additionally take `placement` (:outer or :inner), `offset`, `textoffset`, and `thickness`.
"""
function drawplasmid(chr; outfile = "plasmid.png", offset = 0, colors = Dict(),
                     annotations = [], drawingsize = (1000, 1000),
                     internalcolors = Dict(),
                     defaultgenecolor = RGB(.95, .95, .95),
                     ORI = 1,
                     genethickness = 20, textoffset = 35, figureoffset = (0, 0),
                     genegroups = [],
                     title = "")
    xmax, ymax = drawingdimensions(drawingsize)
    chrlength = length(chr.sequence)
    d = min(xmax, ymax) / 2.5
    drawing = Drawing(xmax, ymax, outfile)
    background("white")
    origin()
    translate(figureoffset...)
    fontsize(28)
    textcentered(title, Point(0, -(ymax / 2.2)))
    circle(O, d, :stroke)
    fontsize(22)
    textcentered(string(chrlength) * " bp")
    fontsize(14)
    genes = @genes(chr, CDS)
    ## draw genes
    for (i, gene) in enumerate(genes)
        loc = locus(gene).position
        a1, a2, amid = geneangles(gene, chrlength, ORI)
        if isempty(annotations)
            genetext = gene.product
        else
            genetext = annotations[i]
        end
        if iscomplement(gene)
            gsave()
                setcolor(get(colors, get(gene, :locus_tag, ""), defaultgenecolor))
                sector(d - genethickness, d, a1, a2, :fillpreserve)
            grestore()
            strokepath()
            textradius = d - textoffset
            if -π/2 < amid < π/2
                text(genetext, Point(textradius * cos(amid), textradius * sin(amid)), halign = :right, valign = :center)
            elseif π/2 <= amid <= 3π/2
                text(genetext, Point(textradius * cos(amid), textradius * sin(amid)), halign = :left, valign = :center)
            end
            ## internal colors
            for (r, c) in get(internalcolors, gene, [])
                a1, a2, amid = geneangles(loc.stop - r.stop, loc.stop - r.start, chrlength, ORI)
                gsave()
                    setcolor(c)
                    sector(d-genethickness, d, a1, a2, :fill)
                grestore()
            end
        else
            gsave()
            setcolor(get(colors, get(gene, :locus_tag, ""), defaultgenecolor))
            sector(d, d + genethickness, a1, a2, :fillpreserve)
            grestore()
            strokepath()
            textradius = d + textoffset
            if -π/2 < amid < π/2
                text(genetext, Point(textradius * cos(amid), textradius * sin(amid)), halign = :left, valign = :center)
            elseif π/2 <= amid <= 3π/2
                text(genetext, Point(textradius * cos(amid), textradius * sin(amid)), halign = :right, valign = :center)
            end
            ## internal colors
            for (r, c) in get(internalcolors, gene, [])
                a1, a2, amid = geneangles(loc.start + r.start, loc.start + r.stop, chrlength, ORI)
                gsave()
                    setcolor(c)
                    sector(d, d+genethickness, a1, a2, :fill)
                grestore()
            end
        end
    end
    ## draw group annotations
    fontsize(18)
    for (i, kv) in enumerate(genegroups)
        k, v = kv
        geneids = k
        annotation = get(v, :text, "")
        offset = get(v, :offset, 0.1d)
        textoffset = get(v, :textoffset, 15)
        thickness = get(v, :thickness, 2)
        placement = get(v, :placement, :outer)
        a1, a2, amid = geneangles(locus(genes[geneids][1]).position.start, locus(genes[geneids][end]).position.stop, chrlength, ORI)
        if placement == :inner
            sector(d-offset, d-offset+thickness, a1, a2, :fill)
            textcurvecentered(annotation, amid, d-offset-textoffset, O)
        elseif placement == :outer
            sector(d+offset, d+offset+thickness, a1, a2, :fill)
            textcurvecentered(annotation, amid, d+offset+textoffset, O)
        end
    end
    finish()
end
