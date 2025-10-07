#
#  ResizeMatrices.jl
#  VertexModel
#
#  Function to adjust the size of system matrices

module ResizeMatrices

# Julia packages
using FromFile
using StaticArrays
using SparseArrays

# Local modules
@from "VertexModelContainers.jl" using VertexModelContainers

function resizeMatrices!(params, matrices, nVertsNew, nEdgesNew, nCellsNew)

    # Assume A and B have already been adjusted 

    # Resize other multidimentional matrices in container 
    matrices.Aᵀ                = spzeros(Int64, nVertsNew, nEdgesNew)
    matrices.Ā                 = spzeros(Int64, nEdgesNew, nVertsNew)
    matrices.Āᵀ                = spzeros(Int64, nVertsNew, nEdgesNew)
    matrices.Bᵀ                = spzeros(Int64, nEdgesNew, nCellsNew)
    matrices.B̄                 = spzeros(Int64, nCellsNew, nEdgesNew)
    matrices.B̄ᵀ                = spzeros(Int64, nEdgesNew, nCellsNew)
    matrices.C                 = spzeros(Int64, nCellsNew, nVertsNew)
    matrices.F                 = spzeros(SVector{2,Float64}, nVertsNew, nCellsNew)
    matrices.edgeMidpointLinks = spzeros(SVector{2,Float64}, nCellsNew, nVertsNew)
    
    # Remove components from stored vectors
    resize!(matrices.cellEdgeCount, nCellsNew)
    resize!(matrices.cellVertexOrders, nCellsNew)
    resize!(matrices.cellEdgeOrders, nCellsNew)
    resize!(matrices.cellPositions, nCellsNew)
    resize!(matrices.cellPerimeters, nCellsNew)
    # resize!(matrices.cellOrientedAreas, nCellsNew)
    # resize!(matrices.cellShapeTensor, nCellsNew)
    resize!(matrices.cellAreas, nCellsNew)
    append!(matrices.cellA₀s, fill(params.A₀, nCellsNew - params.nCells))
    append!(matrices.cellL₀s, fill(params.L₀, nCellsNew - params.nCells))
    resize!(matrices.cellTensions, nCellsNew)
    resize!(matrices.cellPressures, nCellsNew)
    # resize!(matrices.cellPerpAxes, nCellsNew)
    resize!(matrices.cellϵs, nCellsNew)
    resize!(matrices.edgeLengths, nEdgesNew)
    resize!(matrices.edgeTangents, nEdgesNew)
    resize!(matrices.edgeMidpoints, nEdgesNew)
    # resize!(matrices.edgeϵs, nEdgesNew)
    append!(matrices.edgeϵs, fill(MMatrix{3,3,Float64}(zeros(3,3)), nEdgesNew - params.nEdges))
    resize!(matrices.boundaryEdges, nEdgesNew)
    resize!(matrices.timeSinceT1, nEdgesNew)
    resize!(matrices.boundaryVertices, nVertsNew)
    resize!(matrices.vertexAreas, nVertsNew)
    resize!(matrices.totalF,nVertsNew)
    resize!(matrices.externalF,nVertsNew)

    # Update stored number of cells, edges, and vertices
    params.nVerts = nVertsNew
    params.nEdges = nEdgesNew
    params.nCells = nCellsNew

    return nothing

end

export resizeMatrices!

end
