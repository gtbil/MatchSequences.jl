import RawArray
import DelimitedFiles

# https://github.com/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_FmIndex.ipynb

# convert all of the below code from python to julia

struct FMIndex
    # F is the the number of times each character occurs
    # in the original string
    F::Dict{UInt8, UInt64}

    # L is the Burrows-Wheeler transform of the original string
    L::Vector{UInt8}

    # SA is the suffix array, corresponding to L
    SA::Vector{UInt64}

    # T is the tally on the number of times each character occurs
    # along L
    # The amount of space can be reduced substantially by messing with this term
    # T::Dict{UInt8, Vector{UInt64}}
    T::Matrix{UInt64}

    # The last thing you need a map between the original chromosomes
    # and parts of chromosomes, and the suffix array

    # the names of the contigs
    N::Vector{String}

    # positions that map to the suffix array
    C::Vector{UInt64}

    # offsets into the orignial contigs
    O::Vector{UInt64}
end

global const σ = UInt8.(collect("\$ACGT"))
global const σ_to_Col = Dict{UInt8, Int8}(ch => i for (i, ch) in enumerate(σ))
global const SA_frac = 32
global const T_frac = 128

# write function to subset SA to SA_frac proportion and return
function subset_SA(SA::Vector{UInt64})
    str_length = length(SA)
    return SA[1:SA_frac:((str_length ÷ SA_frac) * SA_frac + 1)]
end

# write function to subset SA to SA_frac proportion and return
function subset_T(T::Matrix{UInt64})
    str_length = size(T, 1)

    T_new = Matrix{UInt64}(undef, (str_length - 1) ÷ T_frac + 1, length(σ))

    for ch in σ
        T_new[:,σ_to_Col[ch]] .= T[(1:T_frac:(((str_length - 1) ÷ T_frac) * T_frac + 1)), σ_to_Col[ch]]
    end

    return T_new
end

# pretty print an FMIndex
function Base.display(fm::FMIndex)
    display("σ  : " * String(Char.(σ[2:end])))
    for ch in σ
        println(Char.(ch), " ", Int64.(get(fm.F, ch, zeros(UInt64, length(fm.F)))))
    end
    display("BWT: " * String(Char.(fm.L)))
    println("SA : ", Int64.(fm.SA))
    println("N  : ", fm.N)
    println("C  : ", Int64.(fm.C))
    println("O  : ", Int64.(fm.O))
    for ch in σ
        println(Char.(ch), " ", Int64.(fm.T[:, σ_to_Col[ch]]))
    end
end

# write the FM objects to four separate files
function write_fm(fm::FMIndex, basename::String)
    # first write the alphabet
    RawArray.rawrite(σ, basename * ".a")

    # then the number of each char in the original string
    open(basename * ".F", "w") do io
        DelimitedFiles.writedlm(io, map(x -> get(fm.F, x, UInt64(0)), σ))
    end

    RawArray.rawrite(fm.L, basename * ".L", compress = true)

    RawArray.rawrite(fm.SA, basename * ".SA", compress = true)

    # RawArray.rawrite(reduce(hcat, map(x -> get(fm.T, x, zeros(UInt64, length(fm.L))), σ)), basename * ".T", compress = true)
    RawArray.rawrite(fm.T, basename * ".T", compress = true)

    # finally write the last couple objects, that map the SA to contigs
    open(basename * ".N", "w") do io
        DelimitedFiles.writedlm(io, fm.N)
    end

    open(basename * ".C", "w") do io
        DelimitedFiles.writedlm(io, fm.C)
    end

    open(basename * ".O", "w") do io
        DelimitedFiles.writedlm(io, fm.O)
    end
end

function read_fm(basename::String)
    σ = RawArray.raread( basename * ".a") # read the alphabet
    
    F = Dict{UInt8, UInt64}(σ .=> vec(DelimitedFiles.readdlm(basename * ".F", '\t', UInt64, '\n')))
    L = RawArray.raread(basename * ".L")
    SA = RawArray.raread(basename * ".SA")
    T = RawArray.raread(basename * ".T")

    N = vec(DelimitedFiles.readdlm(basename * ".N", '\t', String, '\n'))
    C = vec(DelimitedFiles.readdlm(basename * ".C", '\t', UInt64, '\n'))
    O = vec(DelimitedFiles.readdlm(basename * ".O", '\t', UInt64, '\n'))

    fm = FMIndex(F, L, SA, T, N, C, O)

    return fm
end

# MatchSequences.read_fm("../data/Coker312_sub.fasta")

# StatProfilerHTML.@profilehtml  MatchSequences.read_fm("../data/Coker312_sub.fasta");

struct FMCheckpoints
    0
end