#
#  OrderAroundCell.jl
#  VertexModel
#
#  Created by Christopher Revell on 01/02/2023.
#
#
# Function to return a *circular* array of vertex labels in the correct clockwise order around a cell, 
# and a *circular* array of edge labels in the correct clockwise order around a cell,
# noting that the first vertex in orderedVertices is on the clockwise end of the first edge in orderedEdges
# ie orderedEdges[1] is the most anticlockwise-ward of the full set of edges and vertices

module OrderAroundCell

using LinearAlgebra
using SparseArrays 
using UnPack
using CircularArrays
using FromFile 

@from "VertexModelContainers.jl" using VertexModelContainers

function orderAroundCell(matrices, i)

    @unpack A, B = matrices

    # Find all edges and vertices for cell i
    cellEdges, edgeOrientations = findnz(B[i,:])
    N = length(cellEdges)

    orderedVerticesAroundEdges = zeros(Int64,N,2)
    # orderedVerticesAroundEdges = MMatrix{N,2,Int64}(zeros(Int64,N,2))
    for j=1:N
        verts, orients = findnz(A[cellEdges[j],:])
        orderedVerticesAroundEdges[j,:] = verts[sortperm(edgeOrientations[j].*orients)]
    end
    
    orderedVerts = CircularArray(ones(Int64,N)) # Ordered list of vertices around cell i in clockwise direction 
    orderedEdges = CircularArray(ones(Int64,N)) # Ordered list of edges around cell i in clockwise direction
    
    for j=1:N
        orderedVerts[j] = orderedVerticesAroundEdges[orderedEdges[j-1],2]
        orderedEdges[j] = findfirst(x->x==orderedVerts[j], @view orderedVerticesAroundEdges[:,1])
    end
    orderedEdges .= cellEdges[orderedEdges[0:N-1]]
    
    return orderedVerts, orderedEdges
end

export orderAroundCell

end