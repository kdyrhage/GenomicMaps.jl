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
        elseif atgene > length(genes)
            println(iobuffer, line)
            continue
        end
        if !stripedgene && isnewgene(line)
            ingenes = true
            stripedgene = false
            atgene += 1
            annotation = replace(replace(string(genes[atgene]), r"^ +" => ""), r"\n +" => "\n")
            line = line[1:end - 2] * ">\n<title>" * annotation * "</title></path>\n"
        elseif occursin(r"<path style=\"fill:none;stroke-width:[^;]+;stroke-linecap:[^;]+;stroke-linejoin:miter;stroke:rgb\([^)]+\);stroke-opacity:1;stroke-miterlimit:10;\" d=\"M \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ \"", line)
            # Second half of spanning genes
            annotation = replace(replace(string(genes[atgene]), r"^ +" => ""), r"\n +" => "\n")
            line = line[1:end - 2] * ">\n<title>" * annotation * "</title></path>\n"
        elseif ingenes && occursin(r"<path", line)
            annotation = replace(replace(string(genes[atgene]), r"^ +" => ""), r"\n +" => "\n")
            line = line[1:end - 2] * ">\n<title>" * annotation * "</title></path>\n"
        end

        if ingenes && occursin(r"<g", line)
            stripedgene = true
        elseif ingenes && occursin(r"</g>", line)
            stripedgene = false
        end

        # leaving the genes section, no more annotations
        if ingenes && occursin(r"d=\"M \S+ \S+ C ", line)
            ingenes = false
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
    # Full, non-spanning
    if occursin(r"<path style=\"fill:none;stroke-width:[^;]+;stroke-linecap:[^;]+;stroke-linejoin:[^;]+;stroke:rgb\([^)]+\);stroke-opacity:1;stroke-miterlimit:10;\" d=\"M \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ \"", line)
        return true
    end
    # Short, non-spanning
    if occursin(r"<path style=\"fill:none;stroke-width:[^;]+;stroke-linecap:[^;]+;stroke-linejoin:[^;]+;stroke:rgb\([^)]+\);stroke-opacity:1;stroke-miterlimit:10;\" d=\"M \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ \"", line)
        c = match(r"d=\"M (\S+) (\S+) L \S+ \S+ L \S+ \S+ L \S+ \S+ L (\S+) (\S+)", line).captures
        if c[1] == c[3] && c[2] == c[4]
           return true
        end
    end
    # First half, complement spanning, long
    if occursin(r"<path style=\"fill:none;stroke-width:[^;]+;stroke-linecap:[^;]+;stroke-linejoin:[^;]+;stroke:rgb\([^)]+\);stroke-opacity:1;stroke-miterlimit:10;\" d=\"M \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ \"", line)
        c = match(r"d=\"M (\S+) \S+ L \S+ \S+ L \S+ \S+ L (\S+) \S+", line).captures
        if parse(Float64, c[1]) > parse(Float64, c[2])
            return true
        end
    end
    # First half complement spanning, short
    if occursin(r"<path style=\"fill:none;stroke-width:[^;]+;stroke-linecap:[^;]+;stroke-linejoin:[^;]+;stroke:rgb\([^)]+\);stroke-opacity:1;stroke-miterlimit:10;\" d=\"M \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ \"", line)
        c = match(r"d=\"M (\S+) (\S+) L \S+ \S+ L (\S+) \S+ L \S+ \S+ L (\S+) (\S+)", line).captures
        if c[1] == c[4] && c[2] != c[5] && c[1] > c[3]
           return true
        end
    end
    # First half, non-compliment
    if occursin(r"<path style=\"fill:none;stroke-width:[^;]+;stroke-linecap:[^;]+;stroke-linejoin:[^;]+;stroke:rgb\([^)]+\);stroke-opacity:1;stroke-miterlimit:10;\" d=\"M \S+ \S+ L \S+ \S+ L \S+ \S+ L \S+ \S+ \"", line)
        c = match(r"d=\"M \S+ (\S+) L \S+ \S+ L \S+ \S+ L \S+ (\S+)", line).captures
        if parse(Float64, c[1]) > parse(Float64, c[2])
            return true
        end
    end
    return false
end
