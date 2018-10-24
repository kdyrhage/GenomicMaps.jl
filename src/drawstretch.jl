function drawgenome(chr; kwargs...)
    p = initialise(chr; kwargs...)
    ### Draw intervals
    drawintervals(p)
    ### Draw genes
    setline(p[:arrowwidth])
    fontsize(p[:genetextsize])
    features = p[:features]
    for g in @genes(chr, :feature in features)
        drawgenetext(p, g, p[:genetextfunction](g))
        genecolours = [get(p[:colourmap], c, p[:defaultcolour]) for c in vcat(p[:colourfunction](g))]
        drawgene(p, g, colours = genecolours)
    end
    sethue("black")
    fontsize(28)
    text(chr.name, Point(p[:xmax] / 2, 5), halign = :center, valign = :top)
    finish()
    return nothing
end


"""
    initialise(chr; kwargs...)

Initialise drawing and return a `Dict` containing properties.
"""
function initialise(chr; drawingsize = "A0landscape",
        outfile = "genomap.svg", nbreaks = 30, textangle = -pi/6,
        numberstyle = 'k', defaultcolour = "gray75", arrowwidth = 12,
        genetextsize = 15, colourfunction = x -> x.feature, colourmap = Dict(),
        features = ["CDS", "rRNA", "tRNA"],
        genetextfunction = x -> get(x, :gene, ""))
    (xmax, ymax) = drawingdimensions(drawingsize)
    drawing = Drawing(xmax, ymax, outfile)
    background("white")
    p = Dict()
    p[:genetextsize] = genetextsize
    p[:chromosome] = chr
    p[:drawing] = drawing
    p[:outfile] = outfile
    p[:interval] = Int(round(length(chr.sequence) / nbreaks, RoundUp))
    p[:textangle] = textangle
    breaks = 1:p[:interval]:length(chr.sequence)
    left = Point[]
    right = Point[]
    @inbounds for i in 1:length(breaks)
        push!(left,  Point(0.02 * xmax, i * ymax / (length(breaks) + 1)))
        push!(right, Point(0.93 * xmax, i * ymax / (length(breaks) + 1)))
    end
    lastpoint = between(left[end], right[end], (length(chr.sequence) - breaks[end]) / (p[:interval] - 1))
    p[:breaks] = breaks
    p[:left] = left
    p[:right] = right
    p[:lastpoint] = lastpoint
    p[:numberstyle] = numberstyle
    p[:offset_y] = ymax / 150
    p[:defaultcolour] = defaultcolour
    p[:arrowwidth] = arrowwidth
    p[:xmax] = xmax
    p[:ymax] = ymax
    p[:genetextoffset] = -p[:arrowwidth] - ymax / 20nbreaks
    p[:features] = features
    p[:colourfunction] = colourfunction
    if !isempty(colourmap)
        p[:colourmap] = colourmap
    else
        p[:colourmap] = generatecolours(p)
    end
    ### genetextfunction can be either a function that, given a Gene, returns
    # the text to plot above the gene, or a Symbol showing which property of the
    # gene to plot
    if genetextfunction isa Symbol
        p[:genetextfunction] = g -> getproperty(g, genetextfunction, "")
    else
        p[:genetextfunction] = genetextfunction
    end
    return p
end


function generatecolours(p)
    chr = p[:chromosome]
    features = p[:features]
    plottedgenes = @genes(chr, :feature in features)
    categories = unique(vcat(p[:colourfunction].(plottedgenes)...))
    colours = Dict()
    cmap = Colors.distinguishable_colors(length(categories))
    for i in eachindex(categories)
        colours[categories[i]] = cmap[i]
    end
    return colours
end


# function initialise(chr::Chromosome; kwargs...)
#     (xmax, ymax) = drawingdimensions(get!(kwargs, :drawingsize, "A0landscape"))
#     drawing = Drawing(xmax, ymax, get!(kwargs, :outfile, "genomap.svg"))
#     background("white")
#     kwargs[:chromosome] = chr
#     get!(kwargs, :nbreaks, 35)
#     get!(kwargs, :textangle, -pi/6)
# end


function formatnumber(p, val::Int)
    if p[:numberstyle] == 'k'
        return @sprintf("%.0ikb", val / 1e3)
    elseif p[:numberstyle] == 'm'
        return @sprintf("%.2fMb", val / 1e6)
    elseif p[:numberstyle] == 'n'
        return @sprintf("%i", val)
    else
        return string(val)
    end
end


function inbreak(p, pos::Real)
    @inbounds for n in 1:length(p[:breaks]) - 1
        if pos >= p[:breaks][n] && pos < p[:breaks][n + 1]
            return n
        end
    end
    return length(p[:breaks])
end


"""
    coordinates(p, pos::Real)

Return a `Point` on the genome map for a position `pos` on the genome that is being plotted.
"""
function coordinates(p, pos::Real)
    b = inbreak(p, pos)
    b < length(p[:breaks]) ?
        between(p[:left][b], p[:right][b], (pos - p[:breaks][b]) / (p[:interval] - 1)) :
        between(p[:left][b], p[:lastpoint], (pos - p[:breaks][b]) / (p[:interval] - 1))
end


function genecoordinates(p, gene)
    return (coordinates(p, gene.locus.position.start), coordinates(p, gene.locus.position.stop))
end


function drawgene(p, gene; colours = [p[:defaultcolour]])
    gsave()
    setline(p[:arrowwidth])
    (p1, p2) = genecoordinates(p, gene)
    if p1.y == p2.y # The gene fits on one line
        if iscomplement(gene)
            striped(genearrow, p, p2, p1, colours)
        else
            striped(genearrow, p, p1, p2, colours)
        end
    else    # The gene does not fit on one line
        position = mean(gene.locus.position)
        genebreak = inbreak(p, position)
        genemidpoint = between(p[:left][genebreak], p[:right][genebreak], (position - p[:breaks][genebreak]) / (p[:interval] - 1))
        breakrange = inbreak(p, gene.locus.position.start):inbreak(p, gene.locus.position.stop)
        if iscomplement(gene)
            striped(partialarrow_stop, p, p[:right][first(breakrange)], p1, colours)
            for r in breakrange[2:end-1]
                striped(partialarrow_mid, p, p[:left][r], p[:right][r], colours)
            end
            striped(partialarrow_start, p, p[:left][last(breakrange)], p2, colours)
        else
            striped(partialarrow_start, p, p[:right][first(breakrange)], p1, colours)
            for r in breakrange[2:end-1]
                striped(partialarrow_mid, p, p[:left][r], p[:right][r], colours)
            end
            striped(partialarrow_stop, p, p[:left][last(breakrange)], p2, colours)
        end
    end
    grestore()
    return nothing
end


function drawgeneproperty(p, gene, property = :gene)
    drawgenetext(p, gene, get(gene, property, ""))
end


function drawgenetext(p, gene, s)
    gsave()
    sethue("black")
    position = mean(gene.locus.position)
    genebreak = inbreak(p, position)
    genemidpoint = between(p[:left][genebreak], p[:right][genebreak], (position - p[:breaks][genebreak]) / (p[:interval] - 1))
    textpoint = genemidpoint + (0, p[:genetextoffset])
    setmatrix([cos(p[:textangle]), sin(p[:textangle]), -sin(p[:textangle]), cos(p[:textangle]), textpoint.x, textpoint.y])
    text(s)
    setmatrix([1, 0, 0, 1, 0, 0])
    grestore()
end


function drawintervals(p)
    gsave()
    setline(3)
    sethue("black")
    @inbounds for n in 1:length(p[:breaks]) - 1
        text(formatnumber(p, ((n - 1) * p[:interval]) + 1),
             p[:left][n] + Point(0, p[:offset_y]), halign = :center, valign = :top)
        text(formatnumber(p, n * p[:interval]),
             p[:right][n] + Point(0, p[:offset_y]), halign = :center, valign = :top)
        line(p[:left][n] + Point(1, 0), p[:right][n] - Point(1, 0), :stroke)
    end
    text(formatnumber(p, ((length(p[:breaks]) - 1) * p[:interval]) + 1), p[:left][end] + Point(0, p[:offset_y]), halign = :center, valign = :top)
    text(formatnumber(p, length(p[:chromosome].sequence)), p[:lastpoint] + Point(0, p[:offset_y]), halign = :center, valign = :top)
    line(p[:left][end] + Point(1, 0), p[:lastpoint] - Point(1, 0), :stroke)
    grestore()
end


function drawingdimensions(drawingsize::Tuple{Real, Real}) drawingsize end
function drawingdimensions(drawingsize::AbstractString)
    if occursin(r"\d+x\d+", drawingsize)
        (xmax, ymax) = parse.(Int, match(r"(\d+)x(\d+)", drawingsize).captures)
        return (xmax, ymax)
    end
    if occursin("landscape", drawingsize)
        psize = replace(drawingsize, "landscape" => "")
        ymax, xmax = Luxor.paper_sizes[psize]
    else
        xmax, ymax = Luxor.paper_sizes[drawingsize]
    end
    return (xmax, ymax)
end
