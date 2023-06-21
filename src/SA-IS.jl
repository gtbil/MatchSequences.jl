# thank you to @mikmoore for this solution!
# https://discourse.julialang.org/t/iterator-for-indexes-of-true-values-in-bitvector-without-allocations/100411/4
struct Findall{I}
	itr::I
end

Base.IteratorSize(::Findall) = Base.SizeUnknown()
Base.eltype(x::Findall) = typeof(firstindex(x.itr))

function Base.iterate(x::Findall, current=firstindex(x.itr))
	next = findnext(x.itr, current)
	isnothing(next) && return nothing
	return (next, next + oneunit(next))
end

# build based on the code here:
# https://zork.net/~st/jottings/sais.html
# http://web.stanford.edu/class/archive/cs/cs166/cs166.1196/lectures/04/Small04.pdf


# return a BitArray indicating if each character is an S-Type character

# return a BitArray with True (1) marked if the the current character is equal 
# to the one immediately before it
# and False (0) otherwise


"""
A non-modifying version of `charBeforeIsSame!`.
"""
function charBeforeIsSame(t::Vector{T})::BitVector where {T<:Unsigned}
    # create an un-initialized bitvector
    issame = BitArray{1}(undef, length(t) + 1)

    charBeforeIsSame!(t, issame)

    return issame
end

# this version should be free of allocations

"""
    charBeforeIsSame!(t, issame)

Helper function. Returns a BitVector that is `length(t) + 1` long.

# Arguments:
* `t::Vector{UInt8}`  : Vector of UInt8 characters.
* `issame::BitVector` : BitVector to store the output in (`undef` is OK)
"""
function charBeforeIsSame!(t::Vector{T}, 
                           issame::BitVector)::BitVector where {T<:Unsigned}

    # set the first value to be false (it cannot be equal to the substring 
    # before it)
    issame[1:1] .= 0

    # use this weird slicing syntax to avoid any allocations
    issame[2:length(t)] .= view(t, 2:length(t)) .== view(t, 1:(length(t)-1))

    # set the empty string to always be not equal to char before it
    issame[end:end] .= 0

    return issame
end

# find the indices of L-Type suffixes
# L-type suffixes are LARGER than the suffix to the right

"""
A non-modifying version of `suffixIsLType!`.
"""
function suffixIsLType(t::Vector{T})::BitVector where {T<:Unsigned}
    # get the sameVec
    issame = BitArray{1}(undef, length(t) + 1)
    isltype = BitArray{1}(undef, length(t) + 1)

    suffixIsLType!(t, issame, isltype)

    return isltype
end

# find the indices of L-Type suffixes
# L-type suffixes are LARGER than the suffix to the right

"""
    suffixIsLType!(t, issame, isltype)

Returns a BitVector that is `length(t) + 1` long.

# Arguments:
* `t::Vector{UInt8}`   : Vector of UInt8 characters.
* `issame::BitVector`  : A BitVector of `length(t) + 1`.
                         It will end up indicating if each character is equal 
                            to the one immediately before it. 
                         `suffixIsLType!` will be called internally.
* `isltype::BitVector` : BitVector to store the output in (`undef` is OK). 
                         Will be modified in place, returned for convenience.
"""
function suffixIsLType!(t::Vector{T}, 
                        issame::BitVector,
                        isltype::BitVector)::BitVector where {T<:Unsigned}
    charBeforeIsSame!(t, issame)

    # we know a character must be L type if it is larger than the character
    # immediately to its right
    # so we can reliably mark these as L type substrings
    @views isltype[1:(length(t)-1)] .= t[1:(length(t)-1)] .> t[2:(length(t))]
    
    # the last character before the empty character must be L type
    isltype[length(t):length(t)] .= 1

    # the empty character is guaranteed to be S type
    isltype[end:end] .= 0

    # now move our way up same_vec, and perform the same ops as we go
    for i in length(t):-1:1
        @views isltype[i:i] .= isltype[i] || (issame[i+1] && (isltype[(i+1)]))
    end

    return isltype
end

"""
A non-modifying version of `suffixIsLMS!`.
"""
function suffixIsLMS(isltype::BitVector)::BitVector
    islms = BitArray{1}(undef, length(isltype))

    suffixIsLMS!(isltype, islms)

    return islms
end

"""
    suffixIsLMS!(isltype, islms)

Given a BitVector indicating if each suffix is a left-most S-type or not.

# Arguments:
* `isltype::BitVector`   : BitVector indicating if each suffix is L-type or not.
* `islms::BitVector`     : BitVector to store the output in.
"""
function suffixIsLMS!(isltype::BitVector,
                      islms::BitVector)::BitVector

    # the first suffix will never be a left-most S-type
    islms[1:1] .= 0

    # a suffix is a left-most S-type if it is an S-type and the character before
    # it is an L-type
    # the syntax below says we will set the output vector to true if the 
    # previous character is L-type (true in the first part before `&.`,
    # and the current character is S-type (broadcasted not with .~).
    @views islms[2:end] .= isltype[1:(end-1)] .& .~ isltype[2:(end)]
    
    return islms
end

"""
    suffixIsLMSEqual(t, lms_vec, offset_a, offset_b)

Given the input sequence and a vector indiciating which characters are L-type,
as well as offsets into each vector, will check to see if two LMS substrings are
equal.

Allocates two ints.

# Arguments:
* `t::Vector{UInt8}`   : Vector of UInt8 characters.
* `islms::BitVector`   : BitVector indicating if each suffix is a left-most 
                           S-type or not.
* `offset_a::UInt64`   : Offset into `s` for the first LMS substring.
* `offset_b::UInt64`   : Offset into `s` for the second LMS substring.
"""
function suffixIsLMSEqual(t::Vector{T},
                          islms::BitVector,
                          offset_a::UInt64,
                          offset_b::UInt64)::Bool where {T<:Unsigned}
    # return false if the one offset equals the end of the string
    if (offset_a == length(t) + 1) | (offset_b == length(t) + 1)
        return false
    end

    # find the next left-most s-type for each offset
    next_a = findnext(islms, offset_a + 1) - 1
    next_b = findnext(islms, offset_b + 1) - 1

    # return if they are not the same length
    if (offset_a - next_a) != (offset_b - next_b)
        return false
    end
    
    # try to see if we can do the check without allocations
    # return allequal(@view(t[offset_a:next_a]) .== @view(t[offset_b:next_b]))
    return all(x -> first(x) == last(x),
               zip(@view(t[offset_a:next_a]), @view(t[offset_b:next_b])))
end

# make the buckets
# it seems like it would be best to do this with a bitmask

"""
A non-modifying version of `findBucketSizes!`.
"""
function findBucketSizes(t::Vector{T}; 
                         alphabet=σ::Vector{T})::Vector{UInt64} where {T<:Unsigned}
    bucketsizes = Vector{UInt64}(undef, length(alphabet))

    findBucketSizes!(t, bucketsizes; alphabet=alphabet)

    return bucketsizes
end

"""
    findBucketSizes!(t, bucketsizes; alphabet=σ)

Given the string, `t`, and one temporary vector, `bucketsizes`, as
well as an optional `alphabet` kwar, will find the size of each bucket. The
primary operation is to computer a temporary bit-mask and immediately consume
it, using `count!`
# Arguments:
* `t::Vector{UInt8}`            : Vector of UInt8 characters.
* `bucketsizes::Vector{UInt64}` : Vector to store the output in.
* `alphabet::Vector{UInt8}`     : Vector of UInt8 characters to use as the 
                                  alphabet.
"""
function findBucketSizes!(t::Vector{T}, 
                          bucketsizes::Vector{UInt64}; 
                          alphabet=σ::Vector{T})::Vector{UInt64} where {T<:Unsigned}

    # make a bitmask for each character, and then sum it
    # this way, the bitmask is consumed as it is produced
    
    for ch in alphabet
        @views count!(==(ch), 
                      bucketsizes[ch:ch],
                      t)
    end

    return bucketsizes
end

"""
A non-modifying version of `findBucketHeads`.
"""
function findBucketHeads(bucketsizes::Vector{UInt64})::Vector{UInt64}
    heads = Vector{UInt64}(undef, length(bucketsizes))
    
    findBucketHeads!(bucketsizes, heads)

    return heads
end

"""
    findBucketHeads!(bucketsizes, heads)

A function that finds where a group of characters will _start_ in the resulting
suffix array. Note position `1` is taken up by the empty string, so this will
need to be trimmed out by the caller.
# Arguments:
* `bucketsizes::Vector{UInt64}` : Vector of UInt64 indicating the size of each
                                bucket.
* `heads::Vector{UInt64}``      : Vector to store the output in.
"""
function findBucketHeads!(bucketsizes::Vector{UInt64},
                          heads::Vector{UInt64})::Vector{UInt64}
    # set the first value to be 0 (will be dealt with later)
    heads[1:1] .= 0

    # replace all the following values with the cumulative sum
    # indicating how many characters are in the string sorted up through
    # this point
    @views cumsum!(heads[2:end], bucketsizes[1:(end-1)])

    # add two to each, since place [1] will always be the empty string
    heads .+= 2

    return heads
end

"""
A non-modifying version of `findBucketTails`.
"""
function findBucketTails(bucketsizes::Vector{UInt64})::Vector{UInt64}
    tails = Vector{UInt64}(undef, length(bucketsizes))
    
    findBucketTails!(bucketsizes, tails)
    
    return tails
end

"""
    findBucketTails!(bucketsizes, tails)

A function that finds where a group of characters will _end_ in the resulting
suffix array. Note position `1` is taken up by the empty string, so this will
need to be trimmed out by the caller.
# Arguments:
* `bucketsizes::Vector{UInt64}` : Vector of UInt64 indicating the size of each
                                bucket.
* `tails::Vector{UInt64}``      : Vector to store the output in.
"""
function findBucketTails!(bucketsizes::Vector{UInt64},
                          tails::Vector{UInt64})::Vector{UInt64}
    # count the number of characters seen up until x point
    @views cumsum!(tails, bucketsizes)

    # add one to each, since place [1] will always be the empty string
    # and we need to account for this
    tails .+= 1

    return tails
end

"""
    guessLMSSort!(t, sa, islms, tails)

Insert the LMS substrings into the suffix array. This is done by iterating
over the `islms` BitVector, and only [yielding](https://discourse.julialang.org/t/iterator-for-indexes-of-true-values-in-bitvector-without-allocations/100411/5)
indexes for those that are true.

# Arguments :
* `t::Vector{UInt8}`     : Vector of UInt8 characters.
* `sa::Vector{UInt64}`   : Vector to store the output in.
* `islms::BitVector`     : BitVector indicating which characters are LMS.
* `tails::Vector{UInt64}`: Vector indicating where each bucket ends. THIS is 
                           also modified!
"""
function guessLMSSort!(t::Vector{T}, 
                       sa::Vector{UInt64},
                       islms::BitVector, 
                       tails::Vector{UInt64})::Vector{UInt64} where {T<:Unsigned}

    # set the empty string to be the first suffix
    sa[1:1] .= length(t) + 1

    # iterate along the true values in islms (no allocations?)
    # until we reach the empty string
    # solution from here
    @views for idx in Findall(islms[1:(end-1)])
        sa[tails[t[idx]]:tails[t[idx]]] .= idx
        tails[t[idx]:t[idx]] .-= 1
    end

    #=
    It may be better here to find the sections of LMS suffixes for each of the 
    buckets (alphabets), and then sort them all at once. This would require a 
    set of temp bitmasks for each part of the alphabet, and then basically map
    the indexes of these set bits to the correct places in sa, along with an
    incrementing counter.
    =#
                       
    return sa
end

"""
    induceSortL!(t, sa, isltype, islms, tails)

# Arguments :
* `t::Vector{UInt8}`     : Vector of UInt8 characters.
* `sa::Vector{UInt64}`   : Vector to store the output in.
* `isltype::BitVector`   : BitVector indicating which characters are L-type.
* `islms::BitVector`     : BitVector indicating which characters are LMS.
* `heads::Vector{UInt64}`: Vector indicating where each bucket starts. THIS is 
                           also modified!
"""
function induceSortL!(t::Vector{T}, 
                      sa::Vector{UInt64},
                      isltype::BitVector,
                      islms::BitVector,
                      bucketsizes::Vector{UInt64}, 
                      heads::Vector{UInt64})::Vector{UInt64} where {T<:Unsigned}
    # reset bucket heads
    findBucketHeads!(bucketsizes, heads)

    # iterate along the true values in islms (no allocations?)
    # until we reach the empty string
    # solution from here
    # @views for idx in Iterators.map(first, Iterators.filter(last, pairs(isltype[1:(end-1)] .& islms[2:end])))
    # for (idx_sa, idx_t) in pairs(sa)
    # TODO : we need to try to improve this....
    for idx_t in sa
        @views if idx_t > 1 && isltype[idx_t-1]
            sa[heads[t[idx_t-1]]:heads[t[idx_t-1]]] .= idx_t - 1
            heads[t[idx_t-1]:t[idx_t-1]] .+= 1
        end
    end
                       
    return sa
end

"""
    induceSortS!(t, sa, isltype, islms, tails)

# Arguments :
* `t::Vector{UInt8}`     : Vector of UInt8 characters.
* `sa::Vector{UInt64}`   : Vector to store the output in.
* `isltype::BitVector`   : BitVector indicating which characters are L-type.
* `islms::BitVector`     : BitVector indicating which characters are LMS.
* `tails::Vector{UInt64}`: Vector indicating where each bucket ends. THIS is 
                           also modified!
"""
function induceSortS!(t::Vector{T}, 
                      sa::Vector{UInt64},
                      isltype::BitVector,
                      islms::BitVector, 
                      bucketsizes::Vector{UInt64},
                      tails::Vector{UInt64})::Vector{UInt64} where {T<:Unsigned}
    # reset bucket tails
    findBucketTails!(bucketsizes, tails)

    for idx_t in Iterators.reverse(sa)
        @views if idx_t > 1 && ~ isltype[idx_t-1]
            sa[tails[t[idx_t-1]]:tails[t[idx_t-1]]] .= idx_t - 1
            tails[t[idx_t-1]:t[idx_t-1]] .-= 1
        end
    end
                       
    return sa
end

"""
Need to deal with these allocations...

Here's the basic strategy (split into multiple functions?)
    * figure out which indices in the SA correspond to LMS suffixes
    * determine which of these (adjacent) LMS suffixes are equal
    * sum along this vector to create a new alphabet
    * report the mapping of these new alphabet characters to the original
      alphabet characters
"""
function makeSuffixArraySummary(t::Vector{T},
                                sa::Vector{UInt64},
                                islms::BitVector,
                                pos_in_t::Vector{UInt64},
                                t_new::Vector{UInt64},
                                tempbool1::BitVector,
                                offset::UInt64,
                                nlms::UInt64) where {T<:Unsigned}
    # iterate over the suffix array
    # and determine if adjacent LMS suffixes are equal
    # make a vector that is the same length as the NUMBER of LMS suffixes
    # this will be strided along the indexes of the suffixarray
    # determine which sa indexes are LMS
    for (idx, v) in enumerate(Findall(islms))
        t_new[(idx + offset):(idx + offset)] .= v
    end

    # the first index will always correspond to the empty string
    # so it will always be a new suffix
    # we are setting this to zero because it will be at the (preceding)
    # 1 - 1 = 0 position in the end
    tempbool1[(1 + offset):(1 + offset)] .= 1

    # instantiate a new variable, a view of the suffix array?
    
    # make the vector showing where in the original string each of the LMS
    # suffixes comes from
    @views pos_in_t[(1 + offset):(offset + nlms)] .= sa[t_new[(1 + offset):(offset + nlms)]]

    # now do to the equality check for adjacent LMS suffixes
    @views for (idx, (off1, off2)) in
        enumerate(zip(pos_in_t[(1 + offset):(offset + nlms - 2)],
                      pos_in_t[(2 + offset):(offset + nlms - 1)]))
        tempbool1[(offset + idx + 1):(offset + idx + 1)] .= ~ suffixIsLMSEqual(t, islms, off1, off2)
    end

    # make a vector to output the summary string
    # this will be confined by numlms
    # for now, lets HARDCODE it as a UInt
    # TODO : parameterize this on the number of LMS suffixes

    # finally, do the actual reindexing
    cumsum!(view(t_new, (1 + offset):(offset + nlms)), 
            view(tempbool1, (1 + offset):(offset + nlms)))

    # need to return:
    # the summary string, which is the string of the reindexed LMS substrings
    # the alphabet size, which is the number of unique LMS substrings
    # the summary suffix offsets, which maps from the summary string
    # to positions in t
    return t_new, pos_in_t
end

"""
A benchmarking function that will eventually be used to call the entire
suffix sorting routine.

Space needed:
* t           -  8 bits × length of input
* sa          - 64 bits × length of input
* isltype     -  1 bit  × length of input
* islms       -  1 bit  × length of input
* lmsstempidx1- 64 bits × length of input
* lmsstempidx2- 64 bits × length of input
* tempbool1   -  1 bit  × length of input
___________________________________________
8 + 64 + 1 + 1 + 64 + 64 + 1  = 203 bits × length of input

* bucketsizes - 64 bits × length of alphabet
* heads       - 64 bits × length of alphabet
* tails       - 64 bits × length of alphabet
* σ_in        -  8 bits × length of alphabet  
++++++++++++++++++++++++++++++++++++++++++++

74 bits × length of input + 200 bits × length of alphabet  

For a standard string of length 2_200_000_000 and 5 letter alphabet, this works
out to 20.35 GB for the string + 125 bytes for the alphabet. The amount needed
for the alphabet will grow as the recursive stack grows.

"""
function benchmark(t::Vector{T}, 
                   sa::Vector{UInt64},
                   isltype::BitVector, 
                   islms::BitVector,
                   bucketsizes::Vector{UInt64}, 
                   heads::Vector{UInt64}, 
                   tails::Vector{UInt64},
                   lmstempidx1::Vector{UInt64},
                   lmstempidx2::Vector{UInt64},
                   lmstempbool1::BitVector;
                   σ_in::Vector{T} = UInt8.(collect(1:4))) where {T<:Unsigned}
    # make sure sa is zeroed out to begin with
    # or else the traversal during induce sort will throw weird errors
    # zero out if this is our first iteration
    sa .= 0

    # allocate a few new variables in the recursive case
    # (can we avoid this?)
    #=
    if length(bucketsizes) != length(σ_in)
        bucketsizes = undef(UInt64, length(σ_in))
    end

    if length(heads) != length(σ_in)
        heads = undef(UInt64, length(σ_in))
    end

    if length(tails) != length(σ_in)
        tails = undef(UInt64, length(σ_in))
    end=#

    # build the type map
    # islms is used as a temporary variable here
    # for tracking if adjacent symbols are equal
    MatchSequences.suffixIsLType!( t , islms, isltype )

    # determine which suffixes are LMS
    # use the type vector to find the left-most S-type suffixes
    MatchSequences.suffixIsLMS!( isltype, islms )

    # save the number of LMS suffixes
    startlms = zero(UInt64)
    numlms = UInt64(sum(islms))

    # get the sizes of each bucket
    MatchSequences.findBucketSizes!(t, bucketsizes; alphabet = σ_in)

    # then find their heads and tails
    # MatchSequences.findBucketHeads!(bucketsizes, heads)
    MatchSequences.findBucketTails!(bucketsizes, tails)

    # insert all the LMS suffixes into approximately the right place
    # in the suffix array
    MatchSequences.guessLMSSort!(t, sa, islms, tails)

    # slot in the L-type suffixes next to the LMS suffixes
    MatchSequences.induceSortL!(t, sa, isltype, islms, bucketsizes, heads)

    # now work backwards on the S-type
    MatchSequences.induceSortS!(t, sa, isltype, islms, bucketsizes, tails)

    # the two args that matter are
    # * lmstempidx2 : t_new
    # * lmstempidx1 : pos_in_t
    MatchSequences.makeSuffixArraySummary(t, sa, islms, 
                                          lmstempidx1, lmstempidx2, 
                                          lmstempbool1,
                                          startlms, numlms)
    # now unwind!

    # make the summary suffix array?
    # do we have this already?

    # accurate LMSSort

    # slot in the L-type suffixes next to the LMS suffixes
    MatchSequences.induceSortL!(t, sa, isltype, islms, bucketsizes, heads)

    # now work backwards on the S-type
    MatchSequences.induceSortS!(t, sa, isltype, islms, bucketsizes, tails)

    return isltype, islms, heads, tails, bucketsizes, sa
end

#=
# σ = UInt8.(collect("\$ACGT"))
# σ = UInt8.(collect("abcdefg"))

# varinfo()
using MatchSequences
import StatsBase
import BenchmarkTools 
import Random

Random.seed!(1234)

alpha = UInt8.(collect("\$ACGT"))
t = UInt8.(StatsBase.sample(alpha, 100_000))
t = UInt8.(StatsBase.sample(alpha, 50))

alpha = UInt8.(collect("\$ACGT"))
t = UInt8.(StatsBase.sample(alpha, 100_000))


alpha = UInt8.(collect("abcdefg"))
t = UInt8.(collect("cabbage"))
t = UInt8.(collect("baabaabac"))
t = UInt8.(collect("baabaabacbaabaabac"))

l_vec = BitArray{1}(undef, length(t) + 1)
lms_vec = BitArray{1}(undef, length(t) + 1)
out_vec = Vector{UInt64}(undef, length(alpha))
heads = Vector{UInt64}(undef, length(alpha))
tails = Vector{UInt64}(undef, length(alpha))

numlms = UInt64(sum(MatchSequences.suffixIsLType(t) |> MatchSequences.suffixIsLMS))

# allocate the islms_subvec
lms_subvec1 =  Vector{UInt64}(undef, length(t) + 1)
lms_subvec2 =  similar(lms_subvec1)

lms_boolvec1 = BitArray{1}(undef, length(t) + 1)

# allocate a vector for the suffix array
sa = Vector{UInt64}(undef, length(t) + 1)
# use zeroes for now

alpha_to_col = Dict{UInt8, UInt8}(ch => i for (i, ch) in enumerate(sort(alpha)))
col = sort(collect(values(alpha_to_col)))

# map to t to t′
t′ = map(x -> alpha_to_col[x], t)

out = MatchSequences.benchmark(t′, sa, l_vec, lms_vec, out_vec, heads, tails,
                               lms_subvec1, lms_subvec2, 
                               lms_boolvec1; σ_in = col);

a,b,c,d,e,f = map(x -> Int.(x), out)

println("_" * join(repeat.(Char.(alpha), e)) * "\n" 
        * join(ifelse.(isone.(a[f]), "L" , "S")) * "\n"
        * join(ifelse.(isone.(b[f]), "^" , " ")) * "\n")




BenchmarkTools.@btime MatchSequences.benchmark($t′, $sa, $l_vec, $lms_vec, $out_vec, $heads, $tails,
                                               $lms_subvec1, $lms_subvec2, 
                                               $lms_boolvec1; σ_in = $col);

@profilehtml  for i in 1:1000
       MatchSequences.benchmark(t′, sa, l_vec, lms_vec, out_vec, heads, tails,
                                             lms_subvec1, lms_subvec2,
                                             lms_boolvec1, lms_boolvec2, numlms; σ_in = col);
end

import SuffixArrays
SuffixArrays.suffixsort(t′)

BenchmarkTools.@btime SuffixArrays.suffixsort($t′);

cabbage  282.400 μs (0 allocations: 0 bytes)
         422.100 μs (10 allocations: 7.03 KiB)
100_000    3.652 ms (0 allocations: 0 bytes) (before recursion!)
100_000    7.352 ms (22 allocations: 4.97 MiB)

=#