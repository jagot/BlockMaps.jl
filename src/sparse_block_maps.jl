const SparseBlockMap{T,B} = BlockMap{T,2,B,<:SparseMatrixCSC{B},<:NTuple{2,<:AbstractUnitRange{Int}}}

@inline SparseBlockMap(::Type{B}, block_sizes::Vararg{AbstractVector{Int}, 2}) where B =
    SparseBlockMap(B, map(blockedrange, block_sizes))

function isblockzero(A::SparseBlockMap, I::Integer, J::Integer)
    colnz = nzrange(A.blocks, J)
    isempty(colnz) || I ∉ rowvals(A.blocks)[colnz]
end
isblockzero(A::SparseBlockMap, I::BlockIndex, J::BlockIndex) =
    isblockzero(A, Int(block(I)), Int(block(J)))

function Base.getindex(A::SparseBlockMap{T}, i::Integer, j::Integer) where T
    I,J = findblockindex.(axes(A), (i,j))
    isblockzero(A,I,J) && return zero(T)
    bi,bj = Int.(block.((I,J)))
    i,j = blockindex.((I,J))
    block_getindex(A.blocks[bi,bj], i, j)
end

function SparseBlockMap(::Type{B}, block_sizes::NTuple{2,<:AbstractUnitRange{Int}}) where B
    blocks = spzeros(B, length.(block_sizes)...)
    BlockMap(blocks, block_sizes)
end

# * Application

function LinearAlgebra.mul!(w::AbstractVector{T}, A::SparseBlockMap{T}, v::AbstractVector{T},
                            α::Number=true, β::Number=false) where T
    m,n = size(A)
    size(w,1) == m && size(v,1) == n ||
        throw(DimensionMismatch("Sizes do not agree"))

    ax = axes(A)
    M,N = length.(ax)

    if iszero(β)
        w .= false # To clear NaNs
    else
        w .*= β
    end

    rows = rowvals(A.blocks)
    vals = nonzeros(A.blocks)

    for J = 1:N
        x = view(v, getindex.(ax, (1,J))[2])
        Threads.@threads for nz in nzrange(A.blocks, J)
            I = rows[nz]
            b = vals[nz]
            y = view(w, getindex.(ax, (I,J))[1])
            mul!(y, b, x, α, true)
        end
    end

    w
end

export SparseBlockMap
