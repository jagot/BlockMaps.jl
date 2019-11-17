module BlockMaps

using LinearAlgebra
using SparseArrays
using BlockArrays
import BlockArrays: AbstractBlockSizes, BlockSizes, blocksizes, global2blockindex,
    getblock, setblock!, cumulsizes, globalrange

using LinearMaps
using LazyArrays
using FillArrays

include("block_maps.jl")
include("sparse_block_maps.jl")
include("linear_maps.jl")

end # module
