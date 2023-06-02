module MatchSequences

import ProfileCanvas
import BenchmarkTools
import StatProfilerHTML

include("bwt.jl")
include("fm.jl")
include("reader.jl")

function main(basename = "./test/test.fasta")
    # get the information for this genome
    genome_info = read_fai(basename * ".fai")
    gap_list = map(i -> i.NAME, genome_info)
    gap_number = zeros(length(gap_list))
    
    # process one chromosome at a time
    contigs = Vector{Vector{UInt8}}()
    contig_names = Vector{String}()
    
    i = 1
    for chr in genome_info
        #println(chr.NAME)
        seq_raw = read_chr(basename, chr)
        pos_Ns = find_Ns(seq_raw)
        gap_number[i] = length(pos_Ns)

        # save the contig names
        append!(contig_names, map(x -> string(chr.NAME, ".", x), 1:(length(pos_Ns) + 1)))
        if gap_number == 0
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

    # make ONE FM index with all the sequences - put '$' between them
    seq = reduce((x, y) -> vcat(x, [UInt8('\$')], y), contigs)

    # make the components of the FMIndex then push it
    sa = Sa(seq)
    bwt = bwtViaSa(seq, sa)
    f = rankBwt(bwt).tots
    t = tallyViaBwt(bwt)
    fm = FMIndex(f, bwt, sa, t)

    # print the FM indexes
    # display.(fms)

    write_fm(fm, basename)

    return fm
end

export main

end
