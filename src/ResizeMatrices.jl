#
#  ResizeMatrices.jl
#  VertexModel
#
#  Created by Christopher Revell on 24/01/2024.
#
#
# Function to adjust the size of system matrices

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
    matrices.F                 = fill(SVector{2, Float64}(zeros(2)), (nVertsNew, nCellsNew))
    matrices.externalF         = fill(SVector{2, Float64}(zeros(2)), nVertsNew)
    matrices.totalF            = fill(SVector{2, Float64}(zeros(2)), nVertsNew)
    matrices.edgeMidpointLinks = fill(SVector{2, Float64}(zeros(2)), (nCellsNew, nVertsNew))
    
    # Remove components from stored vectors
    resize!(matrices.cellEdgeCount, nCellsNew)
    resize!(matrices.cellVertexOrders, nCellsNew)
    resize!(matrices.cellEdgeOrders, nCellsNew)
    resize!(matrices.cellPositions, nCellsNew)
    resize!(matrices.cellPerimeters, nCellsNew)
    resize!(matrices.cellOrientedAreas, nCellsNew)
    resize!(matrices.cellAreas, nCellsNew)
    resize!(matrices.cellTensions, nCellsNew)
    resize!(matrices.cellPressures, nCellsNew)
    resize!(matrices.edgeLengths, nEdgesNew)
    resize!(matrices.edgeTangents, nEdgesNew)
    resize!(matrices.edgeMidpoints, nEdgesNew)
    resize!(matrices.boundaryEdges, nEdgesNew)
    resize!(matrices.timeSinceT1, nEdgesNew)
    resize!(matrices.boundaryVertices, nVertsNew)
    resize!(matrices.vertexAreas, nVertsNew)

    # Update stored number of cells, edges, and vertices
    params.nVerts = nVertsNew
    params.nEdges = nEdgesNew
    params.nCells = nCellsNew

    return nothing

end

export resizeMatrices!

end
