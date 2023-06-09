module MatchSequences

import ProfileCanvas
import BenchmarkTools
import StatProfilerHTML
import Logging
import Dates

include("bwt.jl")
include("fm.jl")
include("reader.jl")

function main(basename = "./test/test.fasta")
    # get the information for this genome
    genome_info = read_fai(basename * ".fai")
    gap_list = map(i -> i.NAME, genome_info)
    gap_number = zeros(Int, length(gap_list))
    
    # process one chromosome at a time
    contigs = Vector{Vector{UInt8}}()
    n = Vector{String}()
    o = Vector{UInt64}()
    
    Logging.@debug "Reading in the genome - " * Dates.format(Dates.now(), "HH:MM:SS.ssss")
    i = 1
    for chr in genome_info
        seq_raw = read_chr(basename, chr)
        pos_Ns = find_Ns(seq_raw)
        gap_number[i] = length(pos_Ns)

        # save the contig names
        append!(n, repeat([chr.NAME],  gap_number[i] + 1))

        push!(o, UInt64(0))

        append!(o, cumsum(map(x -> x.second - x.first + 1, pos_Ns)))

        if gap_number[i] == 0
            push!(contigs, seq_raw)
        else
            for contig in break_on_Ns(seq_raw, pos_Ns)
                # break the chromosomes back into contigs
                push!(contigs, contig)
            end
        end
        # increment the original chromosome counter
        i += 1
    end
 
    # get the map that will also be stored in the FM index
    c = map(x -> UInt64(length(x)), contigs)
    c .+= collect(0:(length(c) - 1))
    cumsum!(c, c)

    # Dates.now make offsets into each chromosome -
    # this is the gap size within each chromosome

    # make ONE FM index with all the sequences - put '$' between them
    seq = reduce((x, y) -> vcat(x, [UInt8('\$')], y), contigs)
    return 0
    Logging.@debug "Making the suffix array - " * Dates.format(Dates.now(), "HH:MM:SS.ssss")
    # make the components of the FMIndex then push it
    sa = Sa(seq)

    Logging.@debug "Doing the bwt - " * Dates.format(Dates.now(), "HH:MM:SS.ssss")
    bwt = bwtViaSa(seq, sa)

    Logging.@debug "Calculating totals of each character - " * Dates.format(Dates.now(), "HH:MM:SS.ssss")
    f = rankBwt(bwt).tots

    Logging.@debug "Making the tally - " * Dates.format(Dates.now(), "HH:MM:SS.ssss")
    t = tallyViaBwt(bwt)

    Logging.@debug "Making the FM index struct - " * Dates.format(Dates.now(), "HH:MM:SS.ssss")

    fm = FMIndex(f, bwt, subset_SA(sa), subset_T(t), n, c, o)
    # print the FM indexes
    # display.(fms)

    Logging.@debug "Writing the FM Index - " *Dates.format(Dates.now(), "HH:MM:SS.ssss")
    write_fm(fm, basename)

    Logging.@debug "DONE - " *Dates.format(Dates.now(), "HH:MM:SS.ssss")
    return fm
end

export main

end
