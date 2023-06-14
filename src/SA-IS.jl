# import LoopVectorization

# global const σ = UInt8.(collect("\$ACGT"))
global const σ = UInt8.(collect("abcdefg"))
global const σ_to_Col = Dict{UInt8, Int8}(ch => i for (i, ch) in enumerate(σ))

# build based on the code here:
# https://zork.net/~st/jottings/sais.html

# return a BitArray indicating if each character is an S-Type character

# return a BitArray with True (1) marked if the the current character is equal 
# to the one immediately before it
# and False (0) otherwise

"""
A non-modifying version of `charBeforeIsSame!`.
"""
function charBeforeIsSame(t::Vector{UInt8})::BitVector
    # create an un-initialized bitvector
    out_vec = BitArray{1}(undef, length(t) + 1)
    return charBeforeIsSame!(t, out_vec)
end

# this version should be free of allocations

"""
    charBeforeIsSame!(t, out_vec)

Helper function. Returns a BitVector that is `length(t) + 1` long.

# Arguments:
* `t::Vector{UInt8}`  : Vector of UInt8 characters.
* `out_vec::BitVector`: BitVector to store the output in (`undef` is OK)
"""
function charBeforeIsSame!(t::Vector{UInt8}, out_vec::BitVector)::BitVector

    # set the first value to be false (it cannot be equal to the substring 
    # before it)
    out_vec[1:1] .= 0

    # use this weird slicing syntax to avoid any allocations
    @views out_vec[2:length(t)] .= t[2:length(t)] .== t[1:(length(t)-1)]

    # set the empty string to always be not equal to char before it
    out_vec[end:end] .= 0

    return out_vec
end

# find the indices of L-Type suffixes
# L-type suffixes are LARGER than the suffix to the right

"""
A non-modifying version of `suffixIsLType!`.
"""
function suffixIsLType(t::Vector{UInt8})::BitVector
    # get the sameVec
    same_vec = BitArray{1}(undef, length(t) + 1)
    out_vec = BitArray{1}(undef, length(t) + 1)

    return suffixIsLType!(t, same_vec, out_vec)
end

# find the indices of L-Type suffixes
# L-type suffixes are LARGER than the suffix to the right

"""
    suffixIsLType!(t, same_vec, out_vec)

Returns a BitVector that is `length(t) + 1` long.

# Arguments:
* `t::Vector{UInt8}`   : Vector of UInt8 characters.
* `same_vec::BitVector`: A BitVector of `length(t) + 1`.
                         It will end up indicating if each character is equal 
                            to the one immediately before it. 
                         `suffixIsLType!` will be called internally.
* `out_vec::BitVector` : BitVector to store the output in (`undef` is OK). 
                         Will be modified in place, returned for convenience.
"""
function suffixIsLType!(t::Vector{UInt8}, 
                        same_vec::BitVector,
                        out_vec::BitVector)::BitVector
    same_vec = charBeforeIsSame!(t, same_vec)

    # we know a character must be L type if it is larger than the character
    # immediately to its right
    # so we can reliably mark these as L type substrings
    @views out_vec[1:(length(t)-1)] .= t[1:(length(t)-1)] .> t[2:(length(t))]
    
    # the last character before the empty character must be L type
    out_vec[length(t):length(t)] .= 1

    # the empty character is guaranteed to be R type
    out_vec[end:end] .= 0

    # now move our way up same_vec, and perform the same ops as we go
    for i in length(t):-1:1
        @views out_vec[i:i] .= out_vec[i] || (same_vec[i+1] && (out_vec[(i+1)]))
    end

    return out_vec
end

"""
A non-modifying version of `suffixIsLMS!`.
"""
function suffixIsLMS(l_vec::BitVector)::BitVector
    out_vec = BitArray{1}(undef, length(l_vec))

    out_vec = suffixIsLMS!(l_vec, out_vec)

    return out_vec
end

"""
Given a BitVector indicating if each suffix is a left-most S-type or not.
# Arguments:
* `l_vec::BitVector`   : BitVector indicating if each suffix is L-type or not.
* `out_vec::BitVector` : BitVector to store the output in.
"""
function suffixIsLMS!(l_vec::BitVector,
                      out_vec::BitVector)::BitVector

    # the first suffix will never be a left-most S-type
    out_vec[1:1] .= 0

    # a suffix is a left-most S-type if it is an S-type and the character before
    # it is an L-type
    # the syntax below says we will set the output vector to true if the 
    # previous character is L-type (true in the first part before `&.`,
    # and the current character is S-type (broadcasted not with .~).
    @views out_vec[2:end] .= l_vec[1:(end-1)] .& .~ l_vec[2:(end)]
    return out_vec
end

"""
    suffixIsLMSEqual(t, lms_vec, offsetA, offsetB)

Given the input sequence and a vector indiciating which characters are L-type,
as well as offsets into each vector, will check to see if two LMS substrings are
equal.

Allocates two ints.

# Arguments:
* `t::Vector{UInt8}`   : Vector of UInt8 characters.
* `lms_vec::BitVector` : BitVector indicating if each suffix is a left-most 
                           S-type or not.
* `offsetA::UInt64`    : Offset into `s` for the first LMS substring.
* `offsetB::UInt64`    : Offset into `s` for the second LMS substring.
"""
function suffixIsLMSEqual(t::Vector{UInt8},
                          lms_vec::BitVector,
                          offsetA::UInt64,
                          offsetB::UInt64)::Bool
    # return false if the one offset equals the end of the string
    if (offsetA == length(t) + 1) | (offsetB == length(t) + 1)
        return false
    end

    # find the next left-most s-type for each offset
    nextA = findnext(lms_vec, offsetA + 1) - 1
    nextB = findnext(lms_vec, offsetB + 1) - 1

    # return if they are not the same length
    if (offsetA - nextA) != (offsetB - nextB)
        return false
    end
    
    # try to see if we can do the check without allocations
    return @views allequal(t[offsetA:nextA] .== t[offsetB:nextB])
end

# make the buckets
# it seems like it would be best to do this with a bitmask

"""
A non-modifying version of `findBucketSizes!`.
"""
function findBucketSizes(t::Vector{UInt8}; alphabet=σ)::Vector{UInt64}
    out_vec = Vector{UInt64}(undef, length(alphabet))

    out_vec = findBucketSizes!(t, out_vec; alphabet=alphabet)

    return out_vec
end

"""
    findBucketSizes!(t, out_vec; alphabet=σ)

Given the string, `t`, and two temporary vectors, `temp_vec` and `out_vec`, as
well as an optional `alphabet` kwar, will find the size of each bucket. The
primary operation is to computer a temporary bit-mask and immediately consume
it, using `count!`
# Arguments:
* `t::Vector{UInt8}`   : Vector of UInt8 characters.
* `out_vec::Vector{UInt64}` : Vector to store the output in.
* `alphabet::Vector{UInt8}` : Vector of UInt8 characters to use as the alphabet.
"""
function findBucketSizes!(t::Vector{UInt8}, 
                          out_vec::Vector{UInt64}; alphabet=σ)::Vector{UInt64}

    # make a bitmask for each character, and then sum it
    # this way, the bitmask is consumed as it is produced
    
    for ch in alphabet
        @views count!(==(ch), 
                      out_vec[σ_to_Col[ch]:σ_to_Col[ch]],
                      t)
    end

    return out_vec
end

"""
A non-modifying version of `findBucketHeads`.
"""
function findBucketHeads(bucket_sz::Vector{UInt64})::Vector{UInt64}
    out_vec = Vector{UInt64}(undef, length(bucket_sz))
    out_vec = findBucketHeads!(bucket_sz, out_vec)
    return out_vec
end

"""
    findBucketHeads!(bucket_sz, out_vec)

A function that finds where a group of characters will _start_ in the resulting
suffix array. Note position `1` is taken up by the empty string, so this will
need to be trimmed out by the caller.
# Arguments:
* `bucket_sz::Vector{UInt64}` : Vector of UInt64 indicating the size of each
                                bucket.
* `out_vec::Vector{UInt64}``  : Vector to store the output in.
"""
function findBucketHeads!(bucket_sz::Vector{UInt64},
                          out_vec::Vector{UInt64})::Vector{UInt64}
    # set the first value to be 0 (will be dealt with later)
    out_vec[1:1] .= 0

    # replace all the following values with the cumulative sum
    # indicating how many characters are in the string sorted up through
    # this point
    @views cumsum!(out_vec[2:end], bucket_sz[1:(end-1)])

    # add two to each, since place [1] will always be the empty string
    out_vec .+= 2

    return out_vec
end

"""
A non-modifying version of `findBucketTails`.
"""
function findBucketTails(bucket_sz::Vector{UInt64})::Vector{UInt64}
    out_vec = Vector{UInt64}(undef, length(bucket_sz))
    out_vec = findBucketTails!(bucket_sz, out_vec)
    return out_vec
end

"""
    findBucketTails!(bucket_sz, out_vec)

A function that finds where a group of characters will _end_ in the resulting
suffix array. Note position `1` is taken up by the empty string, so this will
need to be trimmed out by the caller.
# Arguments:
* `bucket_sz::Vector{UInt64}` : Vector of UInt64 indicating the size of each
                                bucket.
* `out_vec::Vector{UInt64}``  : Vector to store the output in.
"""
function findBucketTails!(bucket_sz::Vector{UInt64},
                          out_vec::Vector{UInt64})::Vector{UInt64}
    # count the number of characters seen up until x point
    @views cumsum!(out_vec, bucket_sz)

    # add one to each, since place [1] will always be the empty string
    # and we need to account for this
    out_vec .+= 1

    return out_vec
end

"""
A benchmarking function that will eventually be used to call the entire
suffix sorting routine.
"""
function benchmark(t::Vector{UInt8}, 
                   l_vec::BitVector, lms_vec::BitVector,
                   out_vec2::Vector{UInt64}, 
                   heads::Vector{UInt64}, 
                   tails::Vector{UInt64};
                   σ_in::Vector{UInt8})
    
    MatchSequences.suffixIsLType!( t , lms_vec, l_vec )
    MatchSequences.suffixIsLMS!( l_vec, lms_vec )
    MatchSequences.findBucketSizes!(t, out_vec2; alphabet = σ_in)
    MatchSequences.findBucketHeads!(out_vec2, heads)
    MatchSequences.findBucketTails!(out_vec2, tails)

    return l_vec, lms_vec, heads, tails
end

#=
# varinfo()
using MatchSequences
import StatsBase
import BenchmarkTools 
import Random

Random.seed!(1234)

alpha = UInt8.(collect("\$ACGT"))
t = UInt8.(StatsBase.sample(alpha, 100_000))

alpha = UInt8.(collect("abcdefg"))
t = UInt8.(collect("cabbage"))

l_vec = BitArray{1}(undef, length(t) + 1)
lms_vec = BitArray{1}(undef, length(t) + 1)
out_vec2 = Vector{UInt64}(undef, length(alpha))
heads = Vector{UInt64}(undef, length(alpha))
tails = Vector{UInt64}(undef, length(alpha))

out = MatchSequences.benchmark(t, l_vec, lms_vec, out_vec2, heads, tails; σ_in = alpha);
map(x -> Int.(x), out)

BenchmarkTools.@btime MatchSequences.benchmark($t, $l_vec, $lms_vec, $out_vec2, $heads, $tails; σ_in = $alpha);

100_000  222.304 ms (0 allocations: 0 bytes)

=#