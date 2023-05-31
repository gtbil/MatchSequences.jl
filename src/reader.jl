include("bwt.jl")

# index your genome so that information is available
struct ContigInfo
    NAME::String
    LENGTH::UInt64
    OFFSET::UInt64
    LINEBASES::UInt64
    LINEWIDTH::UInt64
end

function read_fai(fname::String)::Vector{ContigInfo}
    """
    Reads in a .fai file from samtools faidx and returns the chromosome names and their lengths.
    """
    # read all the lines in the `.fai` file
    lines = readlines(fname)
    # nchr = length(lines)

    # create the Vector of ContigInfo for storing the chromosome information
    lengths = Vector{ContigInfo}()
    for line in readlines(fname)
        this_data = split(line, "\t")
        field1 = this_data[1]
        field2 = parse(UInt64, this_data[2])
        field3 = parse(UInt64, this_data[3])
        field4 = parse(UInt64, this_data[4])
        field5 = parse(UInt64, this_data[5])
        push!(lengths, ContigInfo(field1, field2, field3, field4, field5))
    end

    return lengths
end

function read_chr(basename::String, contig_info::ContigInfo)::Vector{UInt8}
    """
    Reads in a single chromosome from the genome. Return:
        * The chromosome in a Vector{Char}
    Will modify the file-pointer by advancing it.
    """
    # get a file pointer to the `.fasta` file
    fp = open(basename)

    # create a vector of the correct length
    seq = Array{UInt8}(undef, contig_info.LENGTH)
    
    # seek the fp to the OFFSET
    seek(fp, contig_info.OFFSET)
    
    # determine the number of rows we should read "normally" (minus the last one)
    num_lines = contig_info.LENGTH ÷ contig_info.LINEBASES

    # determine the usual line length
    len_usual = contig_info.LINEBASES
    
    # determine the length of the last row
    len_last = contig_info.LENGTH % contig_info.LINEBASES
    
    # determine how many lines in the file to read
    for ptr in 1:num_lines
        seq[((ptr - 1)*len_usual + 1):(ptr*len_usual)] .= convert(Vector{UInt8}, collect(readline(fp)))
    end

    # now deal with the last line
    if len_last != 0
        seq[(num_lines*len_usual + 1):(num_lines*len_usual + len_last)] .= convert(Vector{UInt8}, collect(readline(fp)))
    end
    
    close(fp)
    return seq
end

function find_Ns(seq::Vector{UInt8})::Vector{Pair{UInt64, UInt64}}
    """
    Finds the N gaps and their lengths in each contig. 
    Return a vector of Pairs having the First-Last positions for all N's
    """
    # first find all the N's
    pos_Ns = findall(seq .== UInt8('N'))

    # if pos_Ns is empty, return it (since we will not have any gaps)
    if length(pos_Ns) == 0
        return Vector{Pair{UInt64, UInt64}}()
    end

    # find where the N's are non-contiguously numbered
    range_changes = (pos_Ns[2:(length(pos_Ns))] - pos_Ns[1:(length(pos_Ns) - 1)]) .!= 1

    # find where the ranges start
    starts = vcat(first(pos_Ns), pos_Ns[findall(range_changes) .+ 1])

    # and where they end
    ends = vcat(pos_Ns[findall(range_changes)], last(pos_Ns))

    # now zip the starts and ends together
    pairs = map(Pair, starts, ends)

    return pairs
end

function break_on_Ns(seq::Vector{UInt8}, Ns::Vector{Pair{UInt64, UInt64}})::Vector{Vector{UInt8}}
    """
    Returns a vector of contigs, each of which is a vector of UInt8.
    """

    # pre allocated the vectors that will contains each contig
    # get the total chromosome length
    total_length = length(seq)

    # determine the sizes of each contig
    # since scaffolds will never end with gaps, we can safely exclude that possibility (it would be nonsense!)
    contig_ends = vcat(map(x -> x.first, Ns) .- 1, total_length)
    contig_starts =  vcat(1, map(x -> x.second, Ns) .+ 1)
    contig_lengths = contig_ends - contig_starts
    
    contigs = map((i,j) -> seq[i:j], contig_starts, contig_ends)

    return contigs
end