import RawArray

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
    T::Dict{UInt8, Vector{UInt64}}
end

global const σ = UInt8.(collect("\$ACGT"))

# pretty print an FMIndex
function Base.display(fm::FMIndex)
    display("σ  : " * String(Char.(σ[2:end])))
    display("BWT: " * String(Char.(fm.L)))
    println("SA : ", Int64.(fm.SA))
    for ch in σ
        println(Char.(ch), " ", Int64.(get(fm.T, ch, zeros(UInt64, length(fm.L)))))
    end
end

# write the FM objects to four separate files
function write_fm(fm::FMIndex, basename::String)
    # first write the alphabet
    RawArray.rawrite(σ, basename * ".a")

    # then the number of each char in the original string
    RawArray.rawrite(map(x -> get(fm.F, x, UInt64(0)), σ), basename * ".F")

    RawArray.rawrite(fm.L, basename * ".L")

    RawArray.rawrite(fm.L, basename * ".SA", compress = true)

    RawArray.rawrite(reduce(hcat, map(x -> get(fm.T, x, zeros(UInt64, length(fm.L))), σ)), basename * ".T", compress = true)
    
    #=
    open(basename * ".T", "w") do f
        for ch in σ
            write(f, fm.T[ch])
            write(f, "\n")
        end
    end=#
end

function read_fm(basename::String)
    σ = RawArray.raread( basename * ".a") # read the alphabet
    F = Dict{UInt8, UInt64}(s => v for s in σ, v in RawArray.raread(basename * ".F"))
    L = RawArray.raread(basename * ".L")
    SA = RawArray.raread(basename * ".SA")
    T = Dict{UInt8, Vector{UInt64}}(s => v for s in σ, v in eachcol(RawArray.raread(basename * ".T")))

    fm = FMIndex(F, L, SA, T)

    return fm
end

# MatchSequences.read_fm("../data/Coker312_sub.fasta")

# StatProfilerHTML.@profilehtml  MatchSequences.read_fm("../data/Coker312_sub.fasta");

struct FMCheckpoints
    0
end