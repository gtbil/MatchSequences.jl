import SuffixArrays
# source
# https://www.cs.jhu.edu/~langmea/resources/lecture_notes/10_bwt_and_fm_index_v2.pdf

# probably need to define the alphabet somewhere,
# and allow that to be an argument to these other functions

global const σ = UInt8.(collect("\$ACGT"))
global const σ_to_Col = Dict{UInt8, Int8}(ch => i for (i, ch) in enumerate(σ))

# get type stability for SuffixArrays.suffixsort
sais!(T::Vector{UInt8}, SA::Vector{Int64}, fs::Int64, n::Int64, k::Int64, isbwt::Bool)::Vector{Int64} = SuffixArrays.sais(T::Vector{UInt8}, SA::Vector{Int64}, fs::Int64, n::Int64, k::Int64, isbwt::Bool)

# create a struct to hold the output of rankBwt
struct Rank
    ranks::Vector{UInt64}
    tots::Dict{UInt8, UInt64}
end



# now we want to create a writer for FM index.
# ideally this will have four components:
# file 1: genome.F This is the mapping of the alphabet to the number of times 
#                  it occurs.
# file 2: genome.L This is the Burrows-Wheeler transform of the original string.
# file 3: genome.SA This is the suffix array, corresponding to L.
# file 4: genome.T This is the tally on the number of times each character 
#                  occurs.

@inline function rotations(t)
    """ 
    Return list of rotations of input string t
    """
    tt = vcat(t, [0x24], t, [0x24])
    return map(i -> tt[i:(i+length(t))], 1:(length(t)+1))
end

@inline function bwm(t)
    """
    Return lexicographically sorted list of t’s rotations
    """
    return sort(rotations(t))
end

function bwtViaBwm(t)
    """
    Given T, returns BWT(T) by way of the BWM
    """
    return map(x -> last(x), bwm(t))
end

# for profiling: StatProfilerHTML.@profilehtml bwtViaSa(repeat(convert(Vector{UInt8}, collect("ACCCAGTCCCAGTCA")), 1000))
function bwtViaSa(t::Vector{UInt8})
    """
    Given T, returns BWT(T) by way of the suffix array.
    """
    # and at the beginning for lookup purposes
    tt2 = vcat([0x24], t)

    # make an array to store the results of the suffix array via BWT into
    bw = similar(tt2)

    sa = Sa(t)
    # sa = Sa_new(t)

    # map!(x -> tt[x], bw, suffixArray(t))
    map!(x -> tt2[x], bw, sa)

    # display(String(Char.(bw)))
    return bw
end


function bwtViaSa(t::Vector{UInt8}, sa::Vector{UInt64})
    """
    Given T, returns BWT(T) by way of the suffix array.
    This is the second option, where the suffix array is precomputed.
    """
    # and at the beginning for lookup purposes
    tt2 = vcat([0x24], t)

    # make an array to store the results of the suffix array via BWT into
    bw = similar(tt2)

    map!(x -> tt2[x], bw, sa)

    return bw
end

# this is the old suffix array function
function Sa(t::Vector{UInt8})
    """
    Given T, returns the suffix array.
    """
    # put the '$' at the end for sorting purposes
    tt1 = vcat(t, [0x24])

    # make an array to store the results of the suffix array via BWT into
    sa = zeros(UInt64, length(tt1))

    # map!(x -> tt[x], bw, suffixArray(t))
    sa .= SuffixArrays.suffixsort(tt1)

    # display(String(Char.(bw)))
    return sa
end

# this is the new suffix array function - it might (???) be better: tbd
# need to make sure using Int64 directly doesn't cause any problems
# the source code is very confusing (please use comments if you use weird tricks)
# being too general with the types seems to have hurt performance
# https://github.com/JuliaCollections/SuffixArrays.jl/blob/master/src/sais.jl
function Sa_new(t::Vector{UInt8})
    """
    Given T, returns the suffix array.
    """
    # T = UInt64
    # S = Int64
    n = length(t) + 1
    # U = UInt8
    # I = zeros(S, n)

    # put the '$' at the end for sorting purposes
    tt1 = vcat(t, [0x24])

    # make an array to store the results of the suffix array via BWT into
    sa = zeros(Int64, length(tt1))

    # map!(x -> tt[x], bw, suffixArray(t))
    sa = sais!(tt1, sa, 0, n, 256, false)

    # add one to each value, since this is required
    # this is flexibility built into the SuffixArrays package
    sa .+= 1

    # display(String(Char.(bw)))
    return sa
end



function rankBwt(bw::Vector{UInt8})
    """
    Given BWT string bw, return parallel list of B-ranks. Also
    returns tots: map from character to # times it appears.
    """
    # get the alphabet, which is the unique characters in bw
    tots = Dict(vcat(unique(bw), [0x24]) .=> UInt64(0))

    # set the ranks as zero to begin with
    ranks = zeros(UInt64, length(bw))

    # loop over the burrows wheeler transformed string
    for (i, c) in enumerate(bw)
        # increment the correction position
        tots[c] += 1
        # and add one to our string
        ranks[i] = tots[c]
    end
    return Rank(ranks, tots) 
end

function firstCol(tots::Dict{UInt8, UInt64})
    # https://notebook.community/BenLangmead/comp-genomics-class/notebooks/CG_BWT_Reverse
    """
    Return map from character to the range of rows prefixed by the character.
    """
    tots_keys = sort(collect(keys(tots)))

    output = Dict{UInt8, Pair{UInt64, UInt64}}()
    
    totc = UInt64(1)
    
    for c in tots_keys
        output[c] = Pair(UInt64(totc), UInt64(totc + tots[c] - 1))
        totc += tots[c]
    end

    return output
end

function reverseBwt(bw::Vector{UInt8})
    """
    Make T from BWT(T)
    """
    rank = rankBwt(bw)
    col_F = firstCol(rank.tots)
    rowi = 1

    # create the output vector
    #t = similar(bw)
    #t[end] = 0x24
    t = Vector{UInt8}(undef, length(bw) - 1)

    # probably some way to unroll this loop
    for i in (length(bw)-1):-1:1
        c = bw[rowi]
        t[i] = c
        rowi = col_F[c].first + rank.ranks[rowi] - 1
    end
    return t
end

function tallyViaBwt(bw::Vector{UInt8})
    """
    Make tally table t from BWT string bw
    """
    # create a matrix to store output
    tally = Matrix{UInt64}(undef, length(bw), length(σ))

    for ch in σ
        tally[:,σ_to_Col[ch]] .= bw .== ch
    end

    cumsum!(tally, tally; dims = 1)

    return tally
end