module MatchSequences

import ProfileCanvas
import BenchmarkTools
import StatProfilerHTML

include("bwt.jl")
include("fm.jl")
include("reader.jl")

function main(basename = "../data/GhirsutumCoker_698_v1.0.fasta")
    # get the information for this genome
    genome_info = read_fai(basename * ".fai")
    gap_list = map(i -> i.NAME, genome_info)
    gap_number = zeros(length(gap_list))
    
    # process one chromosome at a time
    contigs = Vector{Vector{UInt8}}()
    contig_names = Vector{String}()

    bwts = Vector{Vector{UInt8}}()

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
            push!(bwts, bwtViaSa(seq_raw))
        else
            for contig in break_on_Ns(seq_raw, pos_Ns)
                # break the chromosomes back into contigs
                push!(contigs, contig)

                # now we can try indexing the genome
                push!(bwts, bwtViaSa(contig))
            end
        end

        # increment the original chromosome counter
        i += 1
    end

end

export main

end
