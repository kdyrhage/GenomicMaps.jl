function getsizes(genomes, regions, chrpadding, regionpadding)
    sizes = []
    for chrs in genomes
        m = -chrpadding
        for chr in chrs
            m += chrpadding - regionpadding
            R = get(regions, chr, [eachindex(chr.sequence)])
            for r in R
                m += length(r) + regionpadding
            end
        end
        push!(sizes, m)
    end
    sizes
end

function getpoint()
    lpoint = rpoint + Point((xmax-rightmargin-leftmargin) * relsize(chrpadding), 0)
end

function formatposition(i)
    if i < 1e3
        string(i)
    elseif i < 1e6
        @sprintf "%.1fkb" i/1e3
    elseif i < 1e9
        @sprintf "%.2fMb" i/1e6
    else
        @sprintf "%.2fGb" i/1e9
    end
end

"""
    compareregions(chrs; kwargs...)

Draw a genomic comparison plot for the genomes in `chrs`.
Some useful supported kwargs are:
* outfile           Path to the output. Supported file formats are PNG, SVG, PDF, and EPS.
* drawingsize       Size of the figure, can be a `Tuple{Int,Int}`, or a `String` such as "1000x1000".
* background        Can be `true`/`false` for white/transparent, or a `Colors.jl`-valid color (default false).
* regions           A `Dict` with `GenomicAnnotations.Record`s as keys, and an `Array`s containing `UnitRange`s with chromosome positions that should be included in the plot. Defaults to `eachindex(chr)`.
* drawpositions     Write the chromosomal positions of each region (default false).
* decorations       A `Dict` with `GenomicAnnotations.Record`s as keys, and `Dict`s as values. The nested `Dict`s have `UnitRange`s as keys (representing chromosome positions) and `NamedTuple`s as values. The `NamedTuple`s can have the fields `color`, `thickness`, and `strand (:+, :-, or :X)`.
* matches           A table (e.g. `Matrix` or `DataFrame`) containing BLAST hits, with columns corresponding to the output from BLAST with the parameter "-outfmt 6".
* matchpadding      Padding in pixels between the chromosome and the bands representing matches.
* matchcolor        Color for matches. Can be a single color or an `AbstractVector` containing separate colors for forwards and reverse matches.
* drawgenes         Bool specifying whether or not to draw CDSs. Can also be a `Function` that takesa `Gene` as input and returns a `Bool`.
* drawpartialgenes  Bool specifying whether or not to draw genes that extend past plotted regions (default true).
* colorfunction     Function that takes a `Gene` as input, and gives a color as output. Defaults to `g -> "lightgrey"`
* textfunction      Function that takes a `Gene` as input, and gives a `String` as output. The returned string is shown above the given gene.
* genetextsize      Font size for gene texts (default 12).
* genethickness     Size of the gene.
* colorlegend       
"""
compareregions(chrs::AbstractVector{C}; kwargs...) where C <: GenomicAnnotations.Record = compareregions([vcat(c...) for c in chrs]; kwargs...)
function compareregions(genomes::AbstractVector;
        outfile = "regions.png",
        drawingsize = (1200, 1000),
        background = false,
        chrpadding = 25,
        regionpadding = 25,
        margins = (10, 10, 50, 10),
        names::AbstractArray = [],
        namesfontsize = 18,
        regions::AbstractDict = Dict(),
        drawpositions::Bool = false,
        centers::AbstractDict = Dict(),
        decorations::AbstractDict = Dict(),
        matchpadding = 0,
        matchcolor = (RGBA(0.5,0,0,0.3), RGBA(0,0,0.5,0.3)),
        colorfunction = g -> "lightgrey",
        genelinecolor = g -> "black",
        genelinethickness = 1,
        textfunction = nothing,
        genetextsize = 12,
        genetextoffset = 15,
        drawgenes = false,
        drawpartialgenes = true,
        colorlegend = nothing,
        colorlegendpos = (0.95, 0.5),
        legendtitle = "",
        legend_title_fontsize = 18,
        legend_tick_fontsize = 15,
        legend_band_height = 1,
        legend_band_padding = 0,
        genethickness = 10,
        matches = [])
    xmax, ymax = drawingdimensions(drawingsize)
    drawing = Drawing(xmax, ymax, outfile)
    leftmargin, rightmargin, topmargin, bottommargin = margins
    genomesizes = getsizes(genomes, regions, chrpadding, regionpadding)
    maxgenomesize = maximum(genomesizes)
    ### Set background color
    if background isa Bool
        background ? Luxor.background(colorant"white") : Luxor.background(colorant"transparent")
    else
        Luxor.background(background)
    end
    if length(names) != length(genomes)
        names = [first(chrs).name for chrs in genomes]
    end
    relsize(s::Number) = s / maxgenomesize
    relsize(r) = length(r) / maxgenomesize
    function getchr(chrname)
        for chrs in genomes
            for chr in chrs
                if chrname == chr.name
                    return chr
                end
            end
        end
    end
    function cumpos(chrs, chrname, pos)
        res = -chrpadding
        for chr in chrs
            res += chrpadding - regionpadding
            R = get(regions, chr, [eachindex(chr.sequence)])
            for region in R
                if chr == chrname && pos in region
                    res += regionpadding + pos - first(region)
                    return res
                else
                    res += regionpadding + length(region)
                end
            end
        end
        res
    end
    getpoint(chrname::AbstractString, pos) = getpoint(getchr(chrname), pos)
    function getpoint(chr, pos)
        i = findfirst(chrs -> chr in chrs, genomes)
        y = topmargin + (i-.5) * (ymax - topmargin - bottommargin) / length(genomes)
        j = findfirst(==(chr), genomes[i])
        R = get(regions, chr, [eachindex(chr.sequence)])
        r = findfirst(r -> pos in r, R)
        isnothing(r) && return nothing
        relpos = cumpos(genomes[i], chr, pos)
        lpoint = between(Point(leftmargin, y), Point(xmax - rightmargin, y), (1-relsize(genomesizes[i])) / 2)
        rpoint = between(Point(leftmargin, y), Point(xmax - rightmargin, y), 1 - ((1-relsize(genomesizes[i])) / 2))
        between(lpoint, rpoint, relpos / genomesizes[i])
    end
    function getpoints(chr, p1, p2)
        i = findfirst(chrs -> chr in chrs, genomes)
        y = topmargin + (i-.5) * (ymax - topmargin - bottommargin) / length(genomes)
        j = findfirst(==(chr), genomes[i])
        p = p1 < p2 ? (p1:p2) : (p1:-1:p2)
        R = get(regions, chr, [eachindex(chr.sequence)])
        R = filter(r -> intersect(r, p) > 1, R)
    end
    ## draw matches (before decorations, so that the matches are on the lowest layer)
    mpoint = Point(0, matchpadding)
    for (i, chrs) in enumerate(genomes[2:end])
        chrsnames = [chr.name for chr in chrs]
        prechrs = genomes[i]
        prechrsnames = [chr.name for chr in prechrs]
        df = filter(r -> r[1] in prechrsnames && r[2] in chrsnames, matches)
        for r in eachrow(df)
            if matchcolor isa Union{Tuple, AbstractVector}
                if r[9] <= r[10]
                    setcolor(matchcolor[1])
                else
                    setcolor(matchcolor[2])
                end
            else
                setcolor(matchcolor)
            end
            points = [getpoint(r[1], r[7]),
                    getpoint(r[1], r[8]),
                    getpoint(r[2], r[10]),
                    getpoint(r[2], r[9])]
            if all(!isnothing, points)
                poly(points .+ [mpoint, mpoint, -mpoint, -mpoint], :fill)
            end
        end
    end
    setcolor("black")
    ## draw genome lines and decorations
    for (i, chrs) in enumerate(genomes)
        y = topmargin + (i-.5) * (ymax - topmargin - bottommargin) / length(genomes)
        genomerelsize = relsize(genomesizes[i])
        lpoint = between(Point(leftmargin, y), Point(xmax - rightmargin, y), (1-genomerelsize) / 2)
        ## write name
        fontsize(namesfontsize)
        text(names[i], Point(10, y), halign = :left, valign = :middle)
        for chr in chrs
            R = get(regions, chr, [eachindex(chr.sequence)])
            dec = get(decorations, chr, [])
            for (i, region) in enumerate(R)
                lpoint = getpoint(chr, first(region))
                rpoint = getpoint(chr, region[end])
                setlinecap("butt")
                line(lpoint, rpoint, :stroke)
                ## Write positions
                gsave(); fontsize(8)
                position_y_offset = 20
                if drawpositions
                    if (rpoint - lpoint).x > 60
                        text(formatposition(first(region)), lpoint + Point(0, position_y_offset), halign = :left)
                        text(formatposition(last(region)), rpoint + Point(0, position_y_offset), halign = :right)
                    elseif (rpoint - lpoint).x > 30
                        text(formatposition(first(region)), lpoint + Point(0, position_y_offset), halign = :left)
                        text(formatposition(last(region)), rpoint + Point(0, position_y_offset), halign = :center)
                    else
                        text(formatposition(first(region)), lpoint + Point(0, position_y_offset), halign = :center)
                        text(formatposition(last(region)), rpoint + Point(0, position_y_offset), halign = :center)
                    end
                end
                grestore()
                ## Draw decorations
                for (dregion, d) in dec
                    dregion = intersect(region, dregion)
                    if length(dregion) > 0
                        gsave()
                            setcolor(get(d, :color, "gray"))
                            thickness = get(d, :thickness, 5)
                            dlpoint = between(lpoint, rpoint, first(dregion) / length(region))
                            drpoint = between(lpoint, rpoint, last(dregion) / length(region))
                            strand = get(d, :strand, :X)
                            if strand == :+
                                box(dlpoint, drpoint - Point(0, thickness), :fill)
                            elseif strand == :-
                                box(dlpoint, drpoint + Point(0, thickness), :fill)
                            else
                                box(dlpoint - Point(0, thickness), drpoint + Point(0, thickness), :fill)
                            end
                        grestore()
                    end
                end
                ## Draw genes if applicable
                fontsize(genetextsize)
                setlinecap(:square)
                if drawgenes != false
                    p = (;arrowwidth = genethickness, genelinethickness = genelinethickness)
                    drawfunction(gene) = drawgenes isa Bool ? feature(gene) == :CDS : drawgenes(gene)
                    for gene in @genes(chr, drawfunction(gene), !isempty(intersect(locus(gene).position, region)))
                        start = getpoint(chr, locus(gene).position.start)
                        stop = getpoint(chr, locus(gene).position.stop)
                        !drawpartialgenes && (isnothing(start) || isnothing(stop)) && @goto out
                        box(lpoint - Point(0, 50), rpoint + Point(0, 50), :clip)
                        if locus(gene).strand == '+'
                            isnothing(start) && (start = lpoint - Point(10, 0))
                            isnothing(stop) && (stop = rpoint + Point(10, 0))
                            genearrow(p, start, stop, colorfunction(gene), genelinecolor(gene))
                        elseif locus(gene).strand == '-'
                            isnothing(start) && (start = lpoint - Point(10, 0))
                            isnothing(stop) && (stop = rpoint + Point(10, 0))
                            genearrow(p, stop, start, colorfunction(gene), genelinecolor(gene))
                        end
                        clipreset()
                        @label out
                    end
                    ## Draw gene text
                    for gene in @genes(chr, drawfunction(gene), !isempty(intersect(locus(gene).position, region)))
                        start = getpoint(chr, locus(gene).position.start)
                        stop = getpoint(chr, locus(gene).position.stop)
                        !drawpartialgenes && (isnothing(start) || isnothing(stop)) && @goto out_text
                        if textfunction isa Function
                            t = textfunction(gene)
                            if !isnothing(t) && !isempty(t)
                                isnothing(start) && (start = lpoint)
                                isnothing(stop) && (stop = rpoint)
                                text(t, midpoint(start, stop) - Point(0, genetextoffset), halign = :center, valign = :baseline)
                            end
                        end
                        @label out_text
                    end
                end
                ## Update the left point for the next region
                lpoint = rpoint + Point((xmax-rightmargin-leftmargin) * relsize(regionpadding), 0)
            end
        end
    end
    legend_width = 30
    if !isnothing(colorlegend)
        startp = Point(xmax * colorlegendpos[1], ymax * colorlegendpos[2])
        fontsize(legend_title_fontsize)
        text(legendtitle, startp - Point(-legend_width/2, legend_band_height * length(colorlegend)) - Point(0, 20), halign = :center)
        fontsize(legend_tick_fontsize)
        for (i, c) in enumerate(colorlegend)
            setcolor(c[1])
            p = startp - Point(0, i*(legend_band_height + legend_band_padding))
            rect(p, legend_width, legend_band_height+1, :fill)
            setcolor(colorant"black")
            !isempty(c[2]) && text(string(c[2]), p + Point(legend_width + 10, 0), valign = :middle)
        end
    end
    finish()
end# 
