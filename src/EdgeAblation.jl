#
#  EdgeAblation.jl
#  VertexModel
#
#
#
# Function that ablates a specified edge j and merges adjacent cells 

module EdgeAblation

using LinearAlgebra
using SparseArrays
using StaticArrays
using DrWatson
using FromFile
using FromFile
using Random
using OrdinaryDiffEq
using InvertedIndices
using Distances

@from "SenseCheck.jl" using SenseCheck
@from "ResizeMatrices.jl" using ResizeMatrices

function edgeAblation!(j, params, matrices, integrator)

    @unpack A,
        B = matrices
    @unpack nVerts,
        nEdges,
        nCells = params 

    # Label adjacent cells as i₁, i₂. We will delete i₂ and add the remaining edges to cell i₁:
    i₁,i₂ = sort(findall(x -> x!=0, B[:,j]))

    cellsToRemove = fill(true,nCells)
    edgesToRemove = fill(true,nEdges)
    vertsToRemove = fill(true,nVerts) # We won't remove vertices in this version of ablation
    
    cellsToRemove[i₂] = false
    edgesToRemove[j] = false

    # Find all edges on i₂, excluding the ablated edge:
    i₂_js = [jj for jj in findall(x -> x!=0, B[i₂,:]) if jj!=j]

    # Add edges to i₁ with same orientation:
    for edge in i₂_js
        B[i₁,edge] = B[i₂,edge] 
    end 

    # Remove edge and cell from adjacency matrices: 
    Btmp = B[cellsToRemove, edgesToRemove]
    Atmp = A[edgesToRemove, vertsToRemove]

    matrices.A = Atmp
    matrices.B = Btmp

    # Resize most matrices in matrices container 
    resizeMatrices!(params, matrices, size(Atmp,2), size(Btmp,2), size(Btmp,1))
    # Some matrices need special treatment because their values cannot be inferred from A, B, and R, so we need to carefully delete specific values
    shrinkIndependentMatrices!(matrices, i₂, findall(x->false, edgesToRemove), [])

    # Update stored number of cells and edges
    params.nEdges = size(Atmp,1)
    params.nCells = size(Btmp,1)


    return nothing

end # end function edgeAblation

function trackVertices!(R,timePoint,params,matrices)

    @unpack trackedVertDistance,
        trackedTimePoints = matrices

    @unpack k_tracked = params

    v1 = R[k_tracked[1]]
    v2 = R[k_tracked[2]]
    dist = sqrt((v1[1]-v2[1])^2 + (v1[2]-v2[2])^2)
    append!(trackedVertDistance,dist)
    append!(trackedTimePoints, timePoint)

    return dist 
    
end # end function trackVertices

export edgeAblation!, trackVertices!

end # end module 

