# thank you to @mikmoore for this solution!
# https://discourse.julialang.org/t/iterator-for-indexes-of-true-values-in-bitvector-without-allocations/100411/4
struct Findall{I}
	itr::I
end

Base.IteratorSize(::Findall) = Base.SizeUnknown()
Base.eltype(x::Findall) = typeof(firstindex(x.itr))

function Base.iterate(x::Findall, current=firstindex(x.itr))
	local next = findnext(x.itr, current)
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

    charBeforeIsSame!(t, issame, UnitRange(UInt64(1), UInt64(length(t)+1)))

    return issame
end

# this version should be free of allocations

"""
charBeforeIsSame!(t, issame)

Helper function. Returns a BitVector that is `length(t) + 1` long.

# Arguments:
* `t::Vector{UInt8}`  : Vector of UInt8 characters.
* `issame::BitVector` : BitVector to store the output in (`undef` is OK)
* `indexes::UnitRange{UInt64}` : Range of indexes to operate on.
"""
function charBeforeIsSame!(t::Vector{T}, 
                           issame::BitVector,
                           indexes::UnitRange{UInt64})::BitVector where {T<:Unsigned}

    # set the first value to be false (it cannot be equal to the substring 
    # before it)
    issame[indexes.start:indexes.start] .= 0

    # use this weird slicing syntax to avoid any allocations
    issame[(indexes.start + 1):(indexes.stop - 1)] .= view(t, (indexes.start + 1):(indexes.stop - 1)) .== view(t, (indexes.start):(indexes.stop - 2))

    # set the empty string to always be not equal to char before it
    issame[indexes.stop:indexes.stop] .= 0

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

    suffixIsLType!(t, issame, isltype, UnitRange(UInt64(1), UInt64(length(t)+1)))

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
                        isltype::BitVector,
                        indexes::UnitRange{UInt64})::BitVector where {T<:Unsigned}
    charBeforeIsSame!(t, issame, indexes)

    # ? local xt = view(t, indexes.start:indexes.stop)
    local xissame = view(issame, indexes)
    local xisltype = view(isltype, indexes)

    # we know a character must be L type if it is larger than the character
    # immediately to its right
    # so we can reliably mark these as L type substrings
    # @views isltype[1:(length(t)-1)] .= t[1:(length(t)-1)] .> t[2:(length(t))]
    @views xisltype[1:(end - 2)] .= t[1:(end - 1)] .> t[2:end]

    # the last character before the empty character must be L type
    xisltype[(end-1):(end-1)] .= 1

    # the empty character is guaranteed to be S type
    xisltype[end:end] .= 0

    # now move our way up same_vec, and perform the same ops as we go
    for i in (length(indexes) - 1):-1:1
        @views xisltype[i:i] .= xisltype[i] || (xissame[i+1] && (xisltype[(i+1)]))
    end

    return isltype
end

"""
A non-modifying version of `suffixIsLMS!`.
"""
function suffixIsLMS(isltype::BitVector)::BitVector
    islms = BitArray{1}(undef, length(isltype))

    suffixIsLMS!(isltype, islms, UnitRange(UInt64(1), UInt64(length(isltype))))

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
                      islms::BitVector,
                      indexes::UnitRange{UInt64})::BitVector
    local xisltype = view(isltype, indexes)
    local xislms = view(islms, indexes)

    # the first suffix will never be a left-most S-type
    xislms[1:1] .= 0

    # a suffix is a left-most S-type if it is an S-type and the character before
    # it is an L-type
    # the syntax below says we will set the output vector to true if the 
    # previous character is L-type (true in the first part before `&.`,
    # and the current character is S-type (broadcasted not with .~).
    @views xislms[2:end] .= xisltype[1:(end-1)] .& .~ xisltype[2:end]
    
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
    local next_a = findnext(islms, offset_a + 1) - 1
    local next_b = findnext(islms, offset_b + 1) - 1

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

    findBucketSizes!(t, bucketsizes,
                     UnitRange(one(UInt64), UInt64(length(t))),
                     UnitRange(one(UInt64), UInt64(length(alphabet))))

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
                          bucketsizes::Vector{UInt64},
                          indexes_t::UnitRange{UInt64},
                          indexes_alpha::UnitRange{UInt64})::Vector{UInt64} where {T<:Unsigned}

    # make a bitmask for each character, and then sum it
    # this way, the bitmask is consumed as it is produced
    local xt = view(t, indexes_t.start:(indexes_t.stop - 1))
    local xbucketsizes = view(bucketsizes, indexes_alpha)
    
    for ch in one(UInt64):length(indexes_alpha)
        count!(==(ch),
               view(xbucketsizes, ch:ch),
               xt)
    end

    return bucketsizes
end


"""
A non-modifying version of `findBucketHeads`.
"""
function findBucketHeads(bucketsizes::Vector{UInt64})::Vector{UInt64}
    heads = Vector{UInt64}(undef, length(bucketsizes))
    
    findBucketHeads!(bucketsizes, heads, 1:length(bucketsizes))

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
                          heads::Vector{UInt64},
                          indexes_alpha::UnitRange{UInt64})::Vector{UInt64}
    local xbucketsizes = view(bucketsizes, indexes_alpha)
    local xheads = view(heads, indexes_alpha)

    # set the first value to be 0 (will be dealt with later)
    xheads[1:1] .= 0

    # replace all the following values with the cumulative sum
    # indicating how many characters are in the string sorted up through
    # this point
    @views cumsum!(xheads[2:end], xbucketsizes[1:(end-1)])

    # add two to each, since place [1] will always be the empty string
    xheads .+= 2

    return heads
end

"""
A non-modifying version of `findBucketTails`.
"""
function findBucketTails(bucketsizes::Vector{UInt64})::Vector{UInt64}
    tails = Vector{UInt64}(undef, length(bucketsizes))
    
    findBucketTails!(bucketsizes, tails, 1:length(bucketsizes))
    
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
                          tails::Vector{UInt64},
                          indexes_alpha::UnitRange{UInt64})::Vector{UInt64}
    local xbucketsizes = view(bucketsizes, indexes_alpha)
    local xtails = view(tails, indexes_alpha)

    # count the number of characters seen up until x point
    cumsum!(xtails, xbucketsizes)

    # add one to each, since place [1] will always be the empty string
    # and we need to account for this
    xtails .+= 1

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
                       tails::Vector{UInt64},
                       indexes_t::UnitRange{UInt64},
                       indexes_alpha::UnitRange{UInt64})::Vector{UInt64} where {T<:Unsigned}
    # make locals for masking
    local xt = view(t, indexes_t.start:(indexes_t.stop - 1))
    local xsa = view(sa, indexes_t)
    local xislms = view(islms, indexes_t)
    local xtails = view(tails, indexes_alpha)

    # set the empty string to be the first suffix
    xsa[1:1] .= length(xt) + 1

    # iterate along the true values in islms (no allocations?)
    # until we reach the empty string
    # solution from here
    @views for idx in Findall(xislms[1:(end-1)])
        xsa[xtails[xt[idx]]:xtails[xt[idx]]] .= idx
        xtails[xt[idx]:xt[idx]] .-= 1
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
                      heads::Vector{UInt64},
                      indexes_t::UnitRange{UInt64},
                      indexes_alpha::UnitRange{UInt64})::Vector{UInt64} where {T<:Unsigned}
    # reset bucket heads
    findBucketHeads!(bucketsizes, heads, indexes_alpha)
    local xsa = view(sa, indexes_t)
    local xheads = view(heads, indexes_alpha)

    # iterate along the true values in islms (no allocations?)
    # until we reach the empty string
    for idx_t in xsa
        @views if idx_t > 1 && isltype[idx_t-1]
            xsa[xheads[t[idx_t-1]]:xheads[t[idx_t-1]]] .= idx_t - 1
            xheads[t[idx_t-1]:t[idx_t-1]] .+= 1
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
                      tails::Vector{UInt64},
                      indexes_t::UnitRange{UInt64},
                      indexes_alpha::UnitRange{UInt64})::Vector{UInt64} where {T<:Unsigned}
    # reset bucket tails
    findBucketTails!(bucketsizes, tails, indexes_alpha)

    local xsa = view(sa, indexes_t)
    local xtails = view(tails, indexes_alpha)

    for idx_t in Iterators.reverse(xsa)
        @views if idx_t > 1 && ~ isltype[idx_t-1]
            xsa[xtails[t[idx_t-1]]:xtails[t[idx_t-1]]] .= idx_t - 1
            xtails[t[idx_t-1]:t[idx_t-1]] .-= 1
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
                                idx_range::UnitRange{UInt64}) where {T<:Unsigned}
    # iterate over the suffix array
    # and determine if adjacent LMS suffixes are equal
    # make a vector that is the same length as the NUMBER of LMS suffixes
    # this will be strided along the indexes of the suffixarray
    # determine which sa indexes are LMS
    for (idx, v) in enumerate(Findall(islms))
        t_new[(idx + idx_range.start - 1):(idx + idx_range.start - 1)] .= v
    end

    # the first index will always correspond to the empty string
    # so it will always be a new suffix
    # we are setting this to zero because it will be at the (preceding)
    # 1 - 1 = 0 position in the end
    tempbool1[idx_range] .= 1

    # instantiate a new variable, a view of the suffix array?
    
    # make the vector showing where in the original string each of the LMS
    # suffixes comes from
    @views pos_in_t[idx_range] .= sa[t_new[idx_range]]

    # now do to the equality check for adjacent LMS suffixes
    @views for (idx, (off1, off2)) in
        enumerate(zip(pos_in_t[idx_range.start:(idx_range.stop - 1)],
                      pos_in_t[(idx_range.start + 1):idx_range.stop]))
        tempbool1[(idx + idx_range.start):(idx + idx_range.start)] .= ~ suffixIsLMSEqual(t, islms, off1, off2)
    end

    # make a vector to output the summary string
    # this will be confined by numlms
    # for now, lets HARDCODE it as a UInt
    # TODO : parameterize this on the number of LMS suffixes

    # finally, do the actual reindexing
    cumsum!(view(t_new, idx_range), 
            view(tempbool1, idx_range))

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
                   t_mapto::Vector{UInt64},
                   t_new::Vector{UInt64},
                   lmstempbool1::BitVector,
                   ranges_sa::Vector{UnitRange{UInt64}},
                   ranges_alpha::Vector{UnitRange{UInt64}};
                   σ_in::Vector{T} = UInt8.(collect(1:4))) where {T<:Unsigned}
    # make sure sa is zeroed out to begin with
    # or else the traversal during induce sort will throw weird errors
    # zero out if this is our first iteration
    sa .= 0
    iter = 1

    # set the initial alphabet size and position in sa
    ranges_alpha[iter] = UnitRange(UInt64(1), UInt64(length(σ_in)))
    ranges_sa[iter] = UnitRange(UInt64(1), UInt64(length(t) + 1))

    # repeat this process until we reach the bottom of the recursion
    # since we doing this as a flat loop, we will track how far along we are
    # with the iter indictor

    # the first pass is different, since it is the only one that operates
    # directly on t

    # islms and isltype will need to keep getting extended
    # as the recursive stack builds

    # length(ranges) is precomputed as ceil(log2(length(t)))
    while iter <= length(ranges_sa)
        # build the type map
        # islms is used as a temporary variable here
        # for tracking if adjacent symbols are equal
        MatchSequences.suffixIsLType!(t, islms, isltype, ranges_sa[iter])

        # determine which suffixes are LMS
        # use the type vector to find the left-most S-type suffixes
        MatchSequences.suffixIsLMS!(isltype, islms, ranges_sa[iter])

        # get the sizes of each bucket
        MatchSequences.findBucketSizes!(t, bucketsizes, ranges_sa[iter], ranges_alpha[iter])

        # then find their heads and tails
        # MatchSequences.findBucketHeads!(bucketsizes, heads)
        MatchSequences.findBucketTails!(bucketsizes, tails, ranges_alpha[iter])

        # insert all the LMS suffixes into approximately the right place
        # in the suffix array
        MatchSequences.guessLMSSort!(t, sa, islms, tails, ranges_sa[iter], ranges_alpha[iter])

        # slot in the L-type suffixes next to the LMS suffixes
        MatchSequences.induceSortL!(t, sa, isltype, islms, bucketsizes, heads, ranges_sa[iter], ranges_alpha[iter])

        # now work backwards on the S-type
        MatchSequences.induceSortS!(t, sa, isltype, islms, bucketsizes, tails, ranges_sa[iter], ranges_alpha[iter])

        #=
        MatchSequences.makeSuffixArraySummary(t, sa, islms, 
                                            t_mapto, t_new, 
                                            lmstempbool1,
                                            ranges_sa[iter])
        # now unwind!
        #println(join(t_mapto, " "))
        #println(join(t_new, " "))

        
        iter += 1
        # accurate LMSSort
        break
    end

    # now do induce sort n times
    while iter >= 1
            # slot in the L-type suffixes next to the LMS suffixes
        MatchSequences.induceSortL!(t, sa, isltype, islms, bucketsizes, heads)

        # now work backwards on the S-type
        MatchSequences.induceSortS!(t, sa, isltype, islms, bucketsizes, tails)
                =#
        break
        
    end

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

l_vec = BitArray{1}(undef, (length(t) + 1)*2)
lms_vec = BitArray{1}(undef, (length(t) + 1)*2)

# this is how deep the "recursive" calls can go in worst case"
ceil_val = ceil(Int64, log2(length(t)))

out_vec      = Vector{UInt64}(undef, length(alpha)^ceil_val)
heads        = Vector{UInt64}(undef, length(alpha)^ceil_val)
tails        = Vector{UInt64}(undef, length(alpha)^ceil_val)
slices_sa    = Vector{UnitRange{UInt64}}(undef, ceil_val)
slices_alpha = Vector{UnitRange{UInt64}}(undef, ceil_val)

numlms = UInt64(sum(MatchSequences.suffixIsLType(t) |> MatchSequences.suffixIsLMS))

# allocate the islms_subvec
lms_subvec1 =  Vector{UInt64}(undef, (length(t) + 1)*2)
lms_subvec2 =  similar(lms_subvec1)

lms_boolvec1 = BitArray{1}(undef, (length(t) + 1)*2)

# allocate a vector for the suffix array
sa = Vector{UInt64}(undef, (length(t) + 1)*2)
# use zeroes for now

alpha_to_col = Dict{UInt8, UInt8}(ch => i for (i, ch) in enumerate(sort(alpha)))
col = sort(collect(values(alpha_to_col)))

# map to t to t′
t′ = map(x -> alpha_to_col[x], t)

out = MatchSequences.benchmark(t′, sa, l_vec, lms_vec, out_vec, heads, tails,
                               lms_subvec1, lms_subvec2, 
                               lms_boolvec1, 
                               slices_sa, slices_alpha; σ_in = col);
#=
a,b,c,d,e,f = map(x -> Int.(x), out)

println("_" * join(repeat.(Char.(alpha), e)) * "\n" 
        * join(ifelse.(isone.(a[f]), "L" , "S")) * "\n"
        * join(ifelse.(isone.(b[f]), "^" , " ")) * "\n")
=#

BenchmarkTools.@btime MatchSequences.benchmark($t′, $sa, $l_vec, $lms_vec, $out_vec, $heads, $tails,
                               $lms_subvec1, $lms_subvec2, 
                               $lms_boolvec1, $slices_sa, $slices_alpha; σ_in = $col);

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