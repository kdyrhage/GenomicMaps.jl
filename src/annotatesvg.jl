function annotatesvg(p)
    ifstream = open(p[:outfile])
    iobuffer = IOBuffer()

    chr = p[:chromosome]
    features = p[:features]
    genes = @genes(chr, feature(gene) in features)

    ingenes = false
    stripedgene = false
    atgene = 0
    atline = 0
    for line in eachline(ifstream)
        atline += 1
        if atline == 2
            line *= "<title>" * p[:chromosome].name * "</title>"
        end
        if isnewgene(line)
            ingenes = true
            stripedgene = false
            atgene += 1
            try
                annotation = replace(replace(string(genes[atgene]), r"^ +" => ""), r"\n +" => "\n")
                line = line[1:end - 2] * ">\n<title>" * annotation * "</title></path>\n"
            catch
            end
        # Second half of spanning genes
    elseif occursin(r"<path style=\"fill-rule:nonzero;fill:rgb(.*);fill-opacity:1;stroke-width:1;stroke-linecap:butt;stroke-linejoin:miter;stroke:rgb(.*);stroke-opacity:1;stroke-miterlimit:10;\" d=\"M \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ \"", line)
            annotation = replace(replace(string(genes[atgene]), r"^ +" => ""), r"\n +" => "\n")
            line = line[1:end - 2] * ">\n<title>" * annotation * "</title></path>\n"
        end

        if ingenes && occursin(r"<g", line)
            stripedgene = true
        elseif ingenes && occursin(r"</g>", line)
            stripedgene = false
        end

        if ingenes && occursin(r"d=\"M \S+ \S+ C ", line)
            ingenes = false
        end

        if stripedgene && occursin(r"<path", line)
            annotation = replace(replace(string(genes[atgene]), r"^ +" => ""), r"\n +" => "\n")
            line = line[1:end - 2] * ">\n<title>" * annotation * "</title></path>\n"
        end

        println(iobuffer, line)
    end

    close(ifstream)

    ofstream = open(p[:outfile], "w")
    print(ofstream, String(take!(iobuffer)))
    close(ofstream)

    return nothing
end


function isnewgene(line)
    if occursin(r"<path style=\"fill-rule:nonzero;fill:rgb(.*);fill-opacity:1;stroke-width:1;stroke-linecap:butt;stroke-linejoin:miter;stroke:rgb(.*);stroke-opacity:1;stroke-miterlimit:10;\" d=\"M \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ \"", line) ||
        # First half of spanning genes on + strand
        occursin(r"<path style=\"fill-rule:nonzero;fill:rgb(.*);fill-opacity:1;stroke-width:1;stroke-linecap:butt;stroke-linejoin:miter;stroke:rgb(.*);stroke-opacity:1;stroke-miterlimit:10;\" d=\"M \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ \"", line)
        return true
    end
        # Second half of spanning genes on - strand
    if occursin(r"<path style=\"fill-rule:nonzero;fill:rgb(.*);fill-opacity:1;stroke-width:1;stroke-linecap:butt;stroke-linejoin:miter;stroke:rgb(.*);stroke-opacity:1;stroke-miterlimit:10;\" d=\"M \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ \"", line)
        c = match(r"d=\"M (\S+) (\S+) L (\S+) (\S+) L (\S+) (\S+) L (\S+) (\S+) L (\S+) (\S+)", line).captures
        if c[2] == c[10]
            return true
        end
    end
    return false
end
