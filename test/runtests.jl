using MatchSequences
using Test

@testset "MatchSequences.jl" begin
    # Write your tests here.
end

#=

### from main
# BenchmarkTools.@btime main("../data/test.fasta")
# StatProfilerHTML.@profilehtml main("../data/test.fasta")
# samtools faidx GhirsutumCoker_698_v1.0.fasta A01 A02 > Coker312_sub.fasta
# samtools faidx Coker312_sub.fasta
# ProfileCanvas.@profview main("../data/Coker312_sub.fasta")
# BenchmarkTools.@btime main("../data/Coker312_sub.fasta")
# main()

# BenchmarkTools.@btime main("../data/test.fasta")
# StatProfilerHTML.@profilehtml main("../data/test.fasta")

# StatProfilerHTML.@profilehtml main("../data/Coker312_sub.fasta")
# BenchmarkTools.@btime main("../data/Coker312_sub.fasta")

import BenchmarkTools
import StatsBase
import Random

Random.seed!(1234)
rand_DNA = convert(Vector{UInt8}, StatsBase.sample(['A', 'C', 'G', 'T'], 20))

# make the SA
SA = Sa(rand_DNA)

# make the bwt
L = bwtViaSa(rand_DNA, SA)

F = rankBwt(L).tots

# use L to create the tallies
T = tallyViaBwt(L)

fm = FMIndex(F, L, SA, T)

display(fm)


# firstCol(rankBwt(bwtViaSa(convert(Vector{UInt8}, collect("abaaba")))).tots)
# BenchmarkTools.@btime Char.(reverseBwt(bwtViaSa(repeat(convert(Vector{UInt8}, collect("abaaba")), 1000000))))
# reverseBwt(bwtViaSa(convert(Vector{UInt8}, collect("ACCCAGTCCCAGTCA"))))
# reverseBwt(bwtViaSa(convert(Vector{UInt8}, collect("abaaba")))

=#