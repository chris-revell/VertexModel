#
#  EdgeAblation.jl
#  VertexModel
#
#  Created by Christopher Revell on 19/12/2023.
#
#
#

module EdgeAblation

# Julia packages
using LinearAlgebra
using SparseArrays
using StaticArrays
using DrWatson
using FromFile
using FromFile
using Random

@from "SenseCheck.jl" using SenseCheck

function edgeAblation(j, params, matrices)
    
    @unpack A, B, cellEdgeCount, cellPositions, cellPerimeters, cellOrientedAreas, cellAreas, cellTensions, cellPressures, cellAges, boundaryEdges, edgeLengths, timeSinceT1, edgeTangents, edgeMidpoints, edgeMidpointLinks, Aᵀ, Ā, Āᵀ, Bᵀ, B̄, B̄ᵀ, C, F, μ, Γ = matrices 
    @unpack nVerts, nEdges, nCells = params

    cells = sort(findall(x->x!=0, B[:,j]))

    matrices.μ[cells[1]] = 0.0
    cell2EdgesTransferredToCell1 = sort([jj for jj in findall(x->x!=0, B[cells[2],:]) if jj!=j])
    cell2EdgesTransferredToCell1Adjusted = Int64[]

    for x=1:length(cell2EdgesTransferredToCell1)
        if cell2EdgesTransferredToCell1[x]>j
            push!(cell2EdgesTransferredToCell1Adjusted,cell2EdgesTransferredToCell1[x]-1)
        else 
            push!(cell2EdgesTransferredToCell1Adjusted,cell2EdgesTransferredToCell1[x])
        end
    end

    newB = B[[ii for ii=1:nCells if ii!=cells[2]],[jj for jj=1:nEdges if jj!=j]]

    newB[cells[1],cell2EdgesTransferredToCell1Adjusted] .= B[cells[2],cell2EdgesTransferredToCell1]

    
    newA = A[[x for x=1:nEdges if x!=j],:]
    
    matrices.A = newA
    matrices.B = newB
            
    deleteat!(cellEdgeCount,cells[2])
    deleteat!(boundaryEdges,cells[2])
    deleteat!(cellPositions,cells[2])
    deleteat!(cellPerimeters,cells[2])
    deleteat!(cellOrientedAreas,cells[2])
    deleteat!(cellAreas,cells[2])
    deleteat!(cellTensions,cells[2])
    deleteat!(cellPressures,cells[2])
    deleteat!(cellAges,cells[2])
    deleteat!(μ,cells[2])
    deleteat!(Γ,j)
    deleteat!(edgeLengths,j)
    deleteat!(edgeTangents,j)
    deleteat!(edgeMidpoints,j)
    deleteat!(timeSinceT1,j)

    matrices.Aᵀ = spzeros(Int64,nVerts,nEdges-1)
    matrices.Ā  = spzeros(Int64,nEdges-1,nVerts)
    matrices.Āᵀ = spzeros(Int64,nVerts,nEdges-1)
    matrices.Bᵀ = spzeros(Int64,nEdges-1,nCells-1)
    matrices.B̄  = spzeros(Int64,nCells-1,nEdges-1)
    matrices.B̄ᵀ = spzeros(Int64,nEdges-1,nCells-1)
    matrices.C  = spzeros(Int64,nCells-1,nVerts)
    matrices.F  = fill(SVector{2,Float64}(zeros(2)), (nVerts,nCells-1))
    matrices.edgeMidpointLinks = fill(SVector{2, Float64}(zeros(2)), (nCells-1, nVerts))

    params.nEdges = nEdges-1
    params.nCells = nCells-1
    
    return nothing 

end

export edgeAblation

end
