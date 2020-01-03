block_eltype(b) = eltype(b)
block_getindex(b, args...) = getindex(b, args...)

struct BlockMap{T, N, B, R<:AbstractArray{B,N}, BS<:NTuple{N,AbstractUnitRange{Int}}} <: AbstractBlockArray{T,N}
    blocks::R
    axes::BS
    BlockMap(blocks::R, block_sizes::BS) where {N, B, R<:AbstractArray{B,N}, BS} =
        new{block_eltype(B), N, B, R, BS}(blocks, block_sizes)
end

@inline Base.axes(block_array::BlockMap) = block_array.axes

isblockzero(A::BlockMap, args...) = false

isblockentryzero(b, args...) = false
isblockentryzero(b::Diagonal, i, j) = i ≠ j
isblockentryzero(b::UpperTriangular, i, j) = i > j
isblockentryzero(b::LowerTriangular, i, j) = i < j

zero_map(::Type{B}, sizes...) where B =
    zeros(eltype(B), sizes...)

@inline function getblock(block_arr::BlockMap{T,N,B}, block::Vararg{Integer, N}) where {T,N,B}
    @boundscheck blockcheckbounds(block_arr, block...)
    if isblockzero(block_arr, block...)
        zero_map(B,blocksize(block_arr.block_sizes, block)...)
    else
        block_arr.blocks[block...]
    end
end

@inline function Base.getindex(block_arr::BlockMap{T,N}, blockindex::BlockIndex{N}) where {T,N}
    @boundscheck blockcheckbounds(block_arr, Block(blockindex.I))
    @inbounds block = getblock(block_arr, blockindex.I...)
    @boundscheck checkbounds(block, blockindex.α...)
    @inbounds v = block[blockindex.α...]
    return v
end


function _check_setblock!(block_arr::BlockMap{T, N}, v, block::NTuple{N, Integer}) where {T,N}
    bls = blocklengths.(axes(block_arr))
    for i in 1:N
        if size(v, i) != bls[i][block[i]]
            throw(DimensionMismatch(string("tried to assign $(size(v)) array to ", getindex.(bls, block), " block")))
        end
    end
end


@inline function setblock!(block_arr::BlockMap{T, N}, v, block::Vararg{Integer, N}) where {T,N}
    @boundscheck blockcheckbounds(block_arr, block...)
    @boundscheck _check_setblock!(block_arr, v, block)
    @inbounds block_arr.blocks[block...] = v
    return block_arr
end

Base.@propagate_inbounds function Base.setindex!(block_array::BlockMap{T, N}, v, block_index::BlockIndex{N}) where {T,N}
    getblock(block_array, block_index.I...)[block_index.α...] = v
end

# * Display

## structured matrix methods ##
function Base.replace_in_print_matrix(A::BlockMap, i::Integer, j::Integer, s::AbstractString)
    I,J = findblockindex.(axes(A), (i,j))
    bi,bj = Int.(block.((I,J)))
    !isblockzero(A,I,J) && !isblockentryzero(A.blocks[bi,bj], blockindex.((I,J))...)  ? s : Base.replace_with_centered_mark(s)
end

block_spy_block(b) = "■"
block_spy_block(::Diagonal) = "⟍"
block_spy_block(::LowerTriangular) = "◣"
block_spy_block(::UpperTriangular) = "◥"

function block_spy(io::IO, A::BlockMap{<:Any,2})
    dm,dn = blocklengths.(axes(A))
    pxl = reduce(gcd, vcat(dm, dn))

    for (I,i) in enumerate(dm)
        for i′ in div(i, pxl)
            for (J,j) in enumerate(dn)
                s = if isblockzero(A, I, J)
                    Base.replace_with_centered_mark(" ")
                else
                    block_spy_block(A.blocks[I,J])
                end
                for j′ in div(j, pxl)
                    write(io, " $s")
                end
            end
            write(io, "\n")
        end
    end    
end

block_spy(A) = block_spy(stdout, A)

export BlockMap
