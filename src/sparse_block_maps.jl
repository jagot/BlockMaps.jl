const SparseBlockMap{T,B} = BlockMap{T,2,B,<:SparseMatrixCSC{B},<:AbstractBlockSizes{2}}

@inline SparseBlockMap(::Type{B}, block_sizes::Vararg{AbstractVector{Int}, 2}) where B =
    SparseBlockMap(B, BlockSizes(block_sizes...))

function isblockzero(A::SparseBlockMap, I, J)
    colnz = nzrange(A.blocks, J)
    isempty(colnz) || I ∉ rowvals(A.blocks)[colnz]
end

function Base.getindex(A::SparseBlockMap{T}, i::Integer, j::Integer) where T
    bi = global2blockindex(A.block_sizes, (i, j))
    I,J = bi.I
    isblockzero(A,I,J) && return zero(T)
    I,J = bi.I
    i,j = bi.α
    block_getindex(A.blocks[I,J], i, j)
end

function SparseBlockMap(::Type{B}, block_sizes::BlockSizes{2}) where B
    n_blocks = nblocks(block_sizes)
    blocks = spzeros(B, n_blocks...)
    BlockMap(blocks, block_sizes)
end

# * Application

function LinearAlgebra.mul!(w::AbstractVector{T}, A::SparseBlockMap{T}, v::AbstractVector{T},
                            α::Number=true, β::Number=false) where T
    m,n = size(A)
    size(w,1) == m && size(v,1) == n ||
        throw(DimensionMismatch("Sizes do not agree"))

    bs = blocksizes(A)
    M,N = nblocks(bs)

    if iszero(β)
        w .= false # To clear NaNs
    else
        w .*= β
    end

    rows = rowvals(A.blocks)
    vals = nonzeros(A.blocks)

    for J = 1:N
        x = view(v, globalrange(bs, (1,J))[2])
        Threads.@threads for nz in nzrange(A.blocks, J)
            I = rows[nz]
            b = vals[nz]
            y = view(w, globalrange(bs, (I,J))[1])
            mul!(y, b, x, α, true)
        end
    end

    w
end

export SparseBlockMap
