const arrowheadangle = 2pi/7
const arrowheadlength = 13
const strokethickness = 1.0


function genearrow(p, start, stop)
    gsave()
    Luxor.translate(start.x, start.y)
    combinedvector = stop - start
    rotate(atan(combinedvector.y, combinedvector.x))
    arrowlength = sqrt(combinedvector.x^2 + combinedvector.y^2)
    points = [
        Point(0, -p[:arrowwidth] / 2),
        Point(arrowlength, 0) + Point(-arrowheadlength * cos(arrowheadangle), -p[:arrowwidth] / 2),
        Point(arrowlength, 0) + Point(-arrowheadlength * cos(arrowheadangle), -arrowheadlength * sin(arrowheadangle)),
        Point(arrowlength, 0),
        Point(arrowlength, 0) + Point(-arrowheadlength * cos(arrowheadangle), arrowheadlength * sin(arrowheadangle)),
        Point(arrowlength, 0) + Point(-arrowheadlength * cos(arrowheadangle), p[:arrowwidth] / 2),
        Point(0, p[:arrowwidth] / 2),
        Point(0, -p[:arrowwidth] / 2)
    ]
    setlinejoin("miter")
    setlinecap("butt")
    poly(points, :fillpreserve)
    sethue("black")
    setline(strokethickness)
    strokepath()
    grestore()
end


function partialarrow_start(p, start, stop)
    gsave()
    Luxor.translate(start.x, start.y)
    combinedvector = stop - start
    rotate(atan(combinedvector.y, combinedvector.x))
    linelength = sqrt(combinedvector.x^2 + combinedvector.y^2)
    points = [
        Point(0, -p[:arrowwidth] / 2),
        Point(linelength, -p[:arrowwidth] / 2),
        Point(linelength, p[:arrowwidth] / 2),
        Point(0, p[:arrowwidth] / 2)
    ]
    poly(points, :fillpreserve)
    sethue("black")
    setline(strokethickness)
    strokepath()
    grestore()
end


function partialarrow_stop(p, start, stop)
    gsave()
    Luxor.translate(start.x, start.y)
    combinedvector = stop - start
    rotate(atan(combinedvector.y, combinedvector.x))
    arrowlength = sqrt(combinedvector.x^2 + combinedvector.y^2)
    points = [
        Point(0, -p[:arrowwidth] / 2),
        Point(arrowlength, 0) + Point(-arrowheadlength * cos(arrowheadangle), -p[:arrowwidth] / 2),
        Point(arrowlength, 0) + Point(-arrowheadlength * cos(arrowheadangle), -arrowheadlength * sin(arrowheadangle)),
        Point(arrowlength, 0),
        Point(arrowlength, 0) + Point(-arrowheadlength * cos(arrowheadangle), arrowheadlength * sin(arrowheadangle)),
        Point(arrowlength, 0) + Point(-arrowheadlength * cos(arrowheadangle), p[:arrowwidth] / 2),
        Point(0, p[:arrowwidth] / 2)
    ]
    setlinejoin("miter")
    setlinecap("butt")
    poly(points, :fillpreserve)
    sethue("black")
    setline(strokethickness)
    strokepath()
    grestore()
end


function partialarrow_mid(p, start, stop)
    gsave()
    Luxor.translate(start.x, start.y)
    combinedvector = stop - start
    rotate(atan(combinedvector.y, combinedvector.x))
    linelength = sqrt(combinedvector.x^2 + combinedvector.y^2)
    points = [Point(0, -p[:arrowwidth] / 2), Point(linelength, -p[:arrowwidth] / 2)]
    poly(points, :stroke)
    strokepath()
    points = [Point(linelength, p[:arrowwidth] / 2), Point(0, p[:arrowwidth] / 2)]
    poly(points, :stroke)
    strokepath()
    grestore()
end


### Striped versions


function stripes(p, start, stop, ncolours, offset = 0)
    stripewidth::Real = 10
    ymax = max(start.y, stop.y) + arrowheadlength
    xmax = max(start.x, stop.x) + arrowheadlength
    function right(point::Real)
        if point <= xmax
            return Point(point, 0)
        else
            return Point(xmax, point - xmax)
        end
    end
    function left(point::Real)
        if point <= xmax
            return Point(0, point)
        else
            return Point(point - ymax, ymax)
        end
    end
    stripepoints = Point[]
    stripeindices = 0:stripewidth:(xmax + ymax + (ncolours * stripewidth))
    currentindex = 1
    while currentindex < div(xmax + ymax, stripewidth)
        append!(stripepoints, (
            right(stripeindices[currentindex]),
            right(stripeindices[currentindex + ncolours]),
            left(stripeindices[currentindex + ncolours]),
            left(stripeindices[currentindex + ncolours + 1])
        ))
        currentindex += ncolours + 1
    end
    append!(stripepoints, (Point(xmax, 0), Point(0, 0)))
    poly(map(point -> point + Point(offset * stripewidth, 0), stripepoints), :clip)
end


function striped(f, p, start, stop, colours; kwargs...)
   gsave()
    for i in 1:length(colours)
        if i > 1
            stripes(p, start, stop, length(colours) - 1, i)
        end
        sethue(colours[i])
        if start.x < stop.x && start.x < (stop.x - (arrowheadlength * cos(arrowheadangle)))
            f(p, start, stop; kwargs...)
        elseif start.x > stop.x && start.x > (stop.x + (arrowheadlength * cos(arrowheadangle)))
            f(p, start, stop; kwargs...)
        elseif start.x < stop.x
            f(p, stop - (arrowheadlength * cos(arrowheadangle), 0), stop; kwargs...)
        elseif start.x > stop.x
            f(p, stop + (arrowheadlength * cos(arrowheadangle), 0), stop; kwargs...)
        end
        clipreset()
    end
    grestore()
end
