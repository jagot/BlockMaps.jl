zero_map(::Type{LinearMap}, M::Integer, N::Integer) =
    LinearMap((y,x) -> (y .= false), M, N,
              ismutating=true, issymmetric=true, ishermitian=true)

block_eltype(::Type{LinearMap{T}}) where T = T
block_spy_block(b::LinearMaps.WrappedMap{<:Any,<:AbstractMatrix}) =
    block_spy_block(b.lmap)

isblockentryzero(b::LinearMaps.WrappedMap{<:Any,<:AbstractMatrix}, args...) =
    isblockentryzero(b.lmap, args...)

unit_vec(::Type{T}, m, i) where T = Vcat(Zeros{T}(i-1), one(T), Zeros{T}(m-i))

function block_getindex(A::LinearMap{T}, i, j) where T
    m,n = size(A)

    ej = unit_vec(T, n, j)
    tmp = Vector{T}(undef, m)

    mul!(tmp, A, ej)

    tmp[i]
end
