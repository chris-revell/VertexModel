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

@from "SenseCheck.jl" using SenseCheck
@from "ResizeMatrices.jl" using ResizeMatrices

function edgeAblation!(j, params, matrices, integrator)

    @unpack A,
        B = matrices
    @unpack nVerts,
        nEdges,
        nCells = params 

    # Find cells adjacent to edge j
    j_is = sort(findall(x -> x!=0, B[:,j]))
    # Find trailing vertices adjacent to edge j left behind after edge ablation 
    j_ks = findall(x -> x!=0, A[j,:])

    cellsToRemove = fill(true,nCells)
    edgesToRemove = fill(true,nEdges)
    vertsToRemove = fill(true,nVerts)
    
    cellsToRemove[j_is] .= false
    edgesToRemove[j] = false

    # Find other edges adjacent to trailing vertices
    for k in j_ks 
       
    end

    return nothing

end # end function edgeAblation

export edgeAblation

end # end module 

