function drawgenome(chr; kwargs...)
    p = initialise(chr; kwargs...)
    ### Draw intervals
    drawintervals(p)
    ### Draw genes
    setline(p[:linewidth])
    fontsize(p[:genetextsize])
    features = p[:features]
    for g in @genes(chr, feature(gene) in features)
        drawgenetext(p, g, p[:genetextfunction](g))
        cv = vcat(p[:colourfunction](g))
        if eltype(cv) <: AbstractFloat
            genecolours = [get(p[:colourmap], c) for c in cv]
        else
            genecolours = [get(p[:colourmap], c, p[:defaultcolour]) for c in cv]
        end
        drawgene(p, g, colours = genecolours)
    end
    ### Draw legend
    drawlegend(p)
    ### Draw title
    sethue("black")
    fontsize(28)
    text(chr.name, Point(p[:xmax] / 2, 5), halign = :center, valign = :top)
    ### Draw extra
    p[:extrafunction](p)
    ### Finish
    finish()
    if get(kwargs, :annotate, true) && occursin(r"\.svg$", p[:outfile])
        annotatesvg(p)
    end
    return p
end


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
    return (coordinates(p, locus(gene).position.start), coordinates(p, locus(gene).position.stop))
end


function drawgene(p, gene; colours = [p[:defaultcolour]])
    gsave()
    setline(p[:genethicknessfunction](gene))
    (p1, p2) = genecoordinates(p, gene)
    if p1.y == p2.y # The gene fits on one line
        if iscomplement(gene)
            striped(genearrow, p, p2, p1, colours)
        else
            striped(genearrow, p, p1, p2, colours)
        end
    else    # The gene does not fit on one line
        position = mean(locus(gene).position)
        genebreak = inbreak(p, position)
        genemidpoint = between(p[:left][genebreak], p[:right][genebreak], (position - p[:breaks][genebreak]) / (p[:interval] - 1))
        breakrange = inbreak(p, locus(gene).position.start):inbreak(p, locus(gene).position.stop)
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
    position = mean(locus(gene).position)
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


function drawlegend(p)
    gsave()
    if p[:legend] == :categorical
        drawlegend_categorical(p)
    elseif p[:legend] == :continuous
        drawlegend_continuous(p)
    else
        return nothing
    end
    grestore()
end

function drawlegend_categorical(p)
    translate(p[:xmax] * p[:rightlimit] * 1.02, p[:ymax] * 0.1)
    setline(3)
    for (category, colour) in p[:colourmap]
        setcolor("black")
        text(string(category), Point(30, 0), valign = :middle)
        box(Point(0, 0), 30, 30, :path)
        strokepreserve()
        setcolor(colour)
        fillpath()
        translate(0, 50)
    end
end


function drawlegend_continuous(p)
    translate(p[:xmax] * p[:rightlimit] * 1.02, p[:ymax] * 0.25)
    setcolor("black")
    fontsize(25)
    text(p[:legend_high], Point(40, 0))
    barheight = p[:ymax] / 2length(p[:colourmap])
    for (i, c) in enumerate(reverse(p[:colourmap]))
        setcolor(c)
        box(Point(0,0), 50, barheight+2, :fill)
        translate(0, barheight)
    end
    setcolor("black")
    text(p[:legend_low], Point(40, 0))
end
