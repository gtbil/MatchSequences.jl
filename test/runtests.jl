import MatchSequences
import Test
import Random
import StatsBase

Random.seed!(1234)

DNA = convert(Vector{UInt8}, StatsBase.sample(['A', 'C', 'G', 'T'], 20))

# run the test fasta
fm_written = MatchSequences.main("./test.fasta")
fm_read = MatchSequences.read_fm("./test.fasta")

Test.@testset "MatchSequences.jl" begin
    # test that Sa and Bw methods for the bwt are the same
    Test.@test MatchSequences.bwtViaBwm(DNA) == MatchSequences.bwtViaSa(DNA)

    # and that they are reversible
    Test.@test MatchSequences.reverseBwt(MatchSequences.bwtViaSa(DNA)) == DNA

    # check which parts of the stucts are equal
    Test.@test fm_written.F == fm_read.F
    Test.@test fm_written.L == fm_read.L
    Test.@test fm_written.SA == fm_read.SA
    Test.@test fm_written.T == fm_read.T
    Test.@test fm_written.N == fm_read.N
    Test.@test fm_written.C == fm_read.C
    Test.@test fm_written.O == fm_read.O
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