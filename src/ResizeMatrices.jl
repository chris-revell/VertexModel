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
    matrices.F                 = spzeros(SVector{2,Float64}, nVertsNew, nCellsNew)
    matrices.edgeMidpointLinks = spzeros(SVector{2,Float64}, nCellsNew, nVertsNew)
    
    # Remove components from stored vectors
    resize!(matrices.cellEdgeCount, nCellsNew)
    resize!(matrices.cellVertexOrders, nCellsNew)
    resize!(matrices.cellEdgeOrders, nCellsNew)
    resize!(matrices.cellPositions, nCellsNew)
    resize!(matrices.cellPerimeters, nCellsNew)
    resize!(matrices.cellOrientedAreas, nCellsNew)
    resize!(matrices.cellShapeTensor, nCellsNew)
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
    resize!(matrices.totalF,nVertsNew)
    resize!(matrices.externalF,nVertsNew)

    resize!(matrices.T_effs,nCellsNew)
    resize!(matrices.P_effs,nCellsNew)
    resize!(matrices.ξs,nCellsNew)
    resize!(matrices.e₁, nCellsNew)
    resize!(matrices.edgeLabels,nEdgesNew)
    resize!(matrices.Λs,nEdgesNew)
    

    # Update stored number of cells, edges, and vertices
    params.nVerts = nVertsNew
    params.nEdges = nEdgesNew
    params.nCells = nCellsNew



    return nothing

end

function shrinkIndependentMatrices!(matrices, removedCells, removedEdges, removedVerts)
    deleteat!(matrices.cellTimeToDivide, removedCells)
    deleteat!(matrices.μ, removedCells)
    deleteat!(matrices.Γ, removedCells)
    deleteat!(matrices.cellA₀s, removedCells)
    deleteat!(matrices.cellL₀s, removedCells)
    deleteat!(matrices.timeSinceT1, removedEdges)
    deleteat!(matrices.firstT1onEdge, removedEdges)
    deleteat!(matrices.cellLabels, removedCells)
    # delete!(matrices.cellsTypeA, removedCells) # Remove the index from the list 
    # delete!(matrices.cellsTypeB, removedCells) # Remove the index from the list 

    
    return nothing 
end

function growIndependentMatrices!(params, matrices, nAddedCells, nAddedEdges)
    append!(matrices.cellTimeToDivide, zeros(nAddedCells))
    append!(matrices.μ, ones(nAddedCells))
    append!(matrices.Γ, fill(params.γ, nAddedCells))
    append!(matrices.cellA₀s, fill(params.A₀, nAddedCells))
    append!(matrices.cellL₀s, fill(params.L₀, nAddedCells))
    append!(matrices.timeSinceT1, fill(params.nonDimCycleTime/100.0, nAddedEdges)) # Ensure new edges can do T1 transitions immediately
    append!(matrices.firstT1onEdge, zeros(nAddedEdges))
    # # In the division case (where a new cell or edge is introduced), these matrices are handled in Division.jl
    # append!(matrices.cellLabels, nCellsNew) 

    return nothing 
end

export resizeMatrices!, shrinkIndependentMatrices!, growIndependentMatrices!

end
