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

    bwts = Vector{Vector{UInt8}}()

    # make a vector for storing the FMIndexes
    fms  = Vector{FMIndex}()
    
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

            # make the components of the FMIndex then push it
            sa = Sa(seq_raw)
            bwt = bwtViaSa(seq_raw, sa)
            f = rankBwt(bwt).tots
            t = tallyViaBwt(bwt)
            push!(fms, FMIndex(f, bwt, sa, t))
            push!(bwts, bwt)
        else
            for contig in break_on_Ns(seq_raw, pos_Ns)
                # break the chromosomes back into contigs
                push!(contigs, contig)

                # make the components of the FMIndex then push it
                sa = Sa(contig)
                bwt = bwtViaSa(contig, sa)
                f = rankBwt(bwt).tots
                t = tallyViaBwt(bwt)
                push!(fms, FMIndex(f, bwt, sa, t))

                # now we can try indexing the genome
                push!(bwts, bwtViaSa(contig))
            end
        end
        # increment the original chromosome counter
        i += 1
    end

    # print the FM indexes
    display.(fms)

    return 0
end

export main

end
