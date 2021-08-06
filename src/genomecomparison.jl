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


"""
    compareregions(chrs; kwargs...)

Draw a genomic comparison plot for the genomes in `chrs`.
Some useful supported kwargs are:
* outfile           Path to the output. Supported file formats are PNG, SVG, PDF, and EPS.
* drawingsize       Size of the figure, can be a `Tuple{Int,Int}`, or a `String` such as "1000x1000".
* regions           A `Dict` with `GenomicAnnotations.Record`s as keys, and an `Array`s containing `UnitRange`s with chromosome positions that should be included in the plot. Defaults to `eachindex(chr)`.
* decorations       A `Dict` with `GenomicAnnotations.Record`s as keys, and `Dict`s as values. The nested `Dict`s have `UnitRange`s as keys (representing chromosome positions) and `NamedTuple`s as values. The `NamedTuple`s can have the fields `color`, `thickness`, and `strand (:+, :-, or :X)`.
* matches           A table (e.g. `Matrix` or `DataFrame`) containing BLAST hits, with columns corresponding to the output from BLAST with the parameter "-outfmt 6".
* matchpadding      Padding in pixels between the chromosome and the bands representing matches.
* drawgenes         Bool specifying whether or not to draw CDSs.
* colorfunction     Function that takes a `Gene` as input, and gives a color as output. Defaults to `g -> "lightgrey"`
* textfunction      Function that takes a `Gene` as input, and gives a `String` as output. The returned string is shown above the given gene.
* genethickness     Size of the gene.
"""
compareregions(chrs::AbstractVector{C}; kwargs...) where C <: GenomicAnnotations.Record = compareregions([vcat(c...) for c in chrs]; kwargs...)
function compareregions(genomes::AbstractVector;
        outfile = "regions.png",
        drawingsize = (1200, 1000),
        chrpadding = 25,
        regionpadding = 25,
        margins = (10, 10, 50, 10),
        names::AbstractArray = [],
        namesfontsize = 18,
        regions::AbstractDict = Dict(),
        centers::AbstractDict = Dict(),
        decorations::AbstractDict = Dict(),
        matchpadding = 0,
        colorfunction = g -> "lightgrey",
        textfunction = nothing,
        genetextoffset = 10,
        drawgenes = false,
        genethickness = 10,
        matches = [])
    xmax, ymax = drawingdimensions(drawingsize)
    drawing = Drawing(xmax, ymax, outfile)
    leftmargin, rightmargin, topmargin, bottommargin = margins
    genomesizes = getsizes(genomes, regions, chrpadding, regionpadding)
    maxgenomesize = maximum(genomesizes)
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
            if r[9] <= r[10]
                setcolor(0.5,0,0,0.3)
            else
                setcolor(0,0,0.5,0.3)
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
                rpoint = lpoint + Point((xmax-rightmargin-leftmargin) * relsize(region), 0)
                line(lpoint, rpoint, :stroke)
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
                                setline(thickness)
                                line(dlpoint, drpoint, :stroke)
                            end
                        grestore()
                    end
                end
                ## Draw genes if applicable
                if drawgenes
                    p = (;arrowwidth = genethickness)
                    for gene in @genes(chr, CDS, !isempty(intersect(locus(gene).position, region)))
                        start = getpoint(chr, locus(gene).position.start)
                        stop = getpoint(chr, locus(gene).position.stop)
                        if locus(gene).strand == '+'
                            if isnothing(start) && isnothing(stop)
                                partialarrow_mid(p, lpoint, rpoint, colorfunction(gene))
                            elseif isnothing(start)
                                partialarrow_stop(p, lpoint, stop, colorfunction(gene))
                            elseif isnothing(stop)
                                partialarrow_start(p, rpoint, start, colorfunction(gene))
                            else
                                genearrow(p, start, stop, colorfunction(gene))
                            end
                        elseif locus(gene).strand == '-'
                            if isnothing(start) && isnothing(stop)
                                partialarrow_mid(p, lpoint, rpoint, colorfunction(gene))
                            elseif isnothing(start)
                                partialarrow_start(p, lpoint, stop, colorfunction(gene))
                            elseif isnothing(stop)
                                partialarrow_stop(p, rpoint, start, colorfunction(gene))
                            else
                                genearrow(p, stop, start, colorfunction(gene))
                            end
                        end
                        ## Draw gene text
                        if textfunction isa Function
                            t = textfunction(gene)
                            if !isnothing(t) && !isempty(t)
                                isnothing(start) && (start = lpoint)
                                isnothing(stop) && (stop = rpoint)
                                text(t, midpoint(start, stop) - Point(0, genetextoffset), halign = :center, valign = :bottom)
                            end
                        end
                    end
                end
                ## Update the left point for the next chr
                lpoint = rpoint + Point((xmax-rightmargin-leftmargin) * relsize(chrpadding), 0)
            end
        end
    end
    finish()
end# 
