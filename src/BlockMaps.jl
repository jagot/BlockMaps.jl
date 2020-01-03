module BlockMaps

using LinearAlgebra
using SparseArrays
using BlockArrays
import BlockArrays: getblock, setblock!, findblockindex, block, blockindex

using LinearMaps
using LazyArrays
using FillArrays

include("block_maps.jl")
include("sparse_block_maps.jl")
include("linear_maps.jl")

export Block

end # module
