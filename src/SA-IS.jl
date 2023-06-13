# build based on the code here:
# https://zork.net/~st/jottings/sais.html

# return a BitArray indicating if each character is an S-Type character


# return a BitArray with True (1) marked if the the current character is equal 
# to the one immediately before it
# and False (0) otherwise

function charBeforeIsSame(t::Vector{UInt8})::BitVector
    # create an un-initialized bitvector
    out_vec = BitArray{1}(undef, length(t) + 1)
    return charBeforeIsSame!(t, out_vec)
end

# this version should be free of allocations
function charBeforeIsSame!(t::Vector{UInt8}, out_vec::BitVector)::BitVector
    # create an un-initialized bitvector

    # set the first value to be false (it cannot be equal to the substring before it)

    # use this weird slicing syntax to avoid any allocations
    out_vec[1:1] .= 0
    @views out_vec[2:length(t)] .= t[2:length(t)] .== t[1:(length(t)-1)]
    out_vec[end:end] .= 0

    return out_vec
end

# find the indices of L-Type suffixes
# L-type suffixes are LARGER than the suffix to the right
function suffixIsLType(t::Vector{UInt8})::BitVector
    # get the sameVec
    same_vec = BitArray{1}(undef, length(t) + 1)
    out_vec = BitArray{1}(undef, length(t) + 1)

    return suffixIsLType!(t, same_vec, out_vec)
end

# find the indices of L-Type suffixes
# L-type suffixes are LARGER than the suffix to the right
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

function suffixIsLMS(l_vec::BitVector)::BitVector
    out_vec = BitArray{1}(undef, length(l_vec))

    out_vec = suffixIsLMS!(l_vec, out_vec)

    return out_vec
end

function suffixIsLMS!(l_vec::BitVector,
                      out_vec::BitVector)::BitVector
    out_vec[1:1] .= 0
    @views out_vec[2:end] .= l_vec[1:(end-1)] .& .~ l_vec[2:(end)]
    return out_vec
end