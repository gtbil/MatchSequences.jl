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

# get the alphabet for an FMIndex
function get_σ(fm::FMIndex)
    return sort(collect(keys(fm.F)))
end

# pretty print an FMIndex
function Base.display(fm::FMIndex)
    σ = get_σ(fm)
    display("σ  : " * String(Char.(σ[2:end])))
    display("BWT: " * String(Char.(fm.L)))
    println("SA : ", Int64.(fm.SA))
    for ch in σ
        println(Char.(ch), " ", Int64.(fm.T[ch]))
    end
end


struct FMCheckpoints
    0
end