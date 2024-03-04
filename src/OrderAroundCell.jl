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


# function orderAroundCell2(matrices, i)

#     @unpack A, B = matrices

#     # Find all edges and vertices for cell i
#     cellEdges, edgeOrientations = findnz(B[i,:])
#     # cellEdges, edgeOrientations = findnz(@view B[i,:])
#     N = length(cellEdges)


#     # a, b, c = findnz(A[cellEdges,:])
#     a, b, c = findnz(@view A[cellEdges,:])
#     b .= b[sortperm(a)]
#     c .= c[sortperm(a)]


#     orderedVerticesAroundEdges = zeros(Int64,N,2)
#     # orderedVerticesAroundEdges = MMatrix{N,2,Int64}(zeros(Int64,N,2))
#     for j=1:N
#         verts = @view b[2*(j-1)+1:2*j]
#         orients = @view c[2*(j-1)+1:2*j]
#         # verts, orients = findnz(A[cellEdges[j],:])
#         orderedVerticesAroundEdges[j,:] .= verts[sortperm(edgeOrientations[j].*orients)]
#     end
    
#     orderedVerts = CircularArray(ones(Int64,N)) # Ordered list of vertices around cell i in clockwise direction 
#     orderedEdges = CircularArray(ones(Int64,N)) # Ordered list of edges around cell i in clockwise direction
    
#     for j=1:N
#         orderedVerts[j] = orderedVerticesAroundEdges[orderedEdges[j-1],2]
#         orderedEdges[j] = findfirst(x->x==orderedVerts[j], @view orderedVerticesAroundEdges[:,1])
#     end
#     orderedEdges .= cellEdges[orderedEdges[0:N-1]]
    
#     return orderedVerts, orderedEdges
# end


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


# function orderAroundCell1(matrices, i)

#     @unpack A, B, C = matrices

#     # Find all edges and vertices for cell i
#     cellEdges = findnz(B[i,:])[1]
#     cellVerts = findnz(C[i,:])[1]
#     N = length(cellEdges)

#     orderedVertices = zeros(Int64,N) # Ordered list of vertices around cell i in clockwise direction 
#     orderedEdges = zeros(Int64,N)    # Ordered list of edges around cell i in clockwise direction
    
#     # Pick a random vertex to start ordering from 
#     nextVertex = [cellVerts[1]]
    
#     for jj = 1:N
#         allNeighbouringEdges = findall(j->j!=0,A[:,nextVertex[end]]) # Find all neighbouring edges around vertex (could be up to 3)                
#         iNeighbourEdges = [j for j in allNeighbouringEdges if B[i,j]!=0] # Find the intersection of all neighbours with edges of cell i to give 2 relevant edges
#         # Testing both edges in iNeighbourEdges
#         # Downstream clockwise if B[i,j]>0 and A[j,k]<0 or B[i,j]<0 and A[j,k]>0
#         if B[i,iNeighbourEdges[1]] > 0 && A[iNeighbourEdges[1],nextVertex[end]] < 0 
#             # iNeighbourEdges[1] is downstream clockwise
#             push!(orderedEdges, iNeighbourEdges[1])
#             push!(orderedVertices, findall(k->k>0,A[iNeighbourEdges[1],:])[1]) # Since A[iNeighbourEdges[1],nextVertex[end]]<0, the next vertex downstream must have A[iNeighbourEdges[1],k]>0
#         elseif B[i,iNeighbourEdges[1]] < 0 && A[iNeighbourEdges[1],nextVertex[end]] > 0 
#             # iNeighbourEdges[1] is downstream clockwise
#             push!(orderedEdges, iNeighbourEdges[1])
#             push!(orderedVertices, findall(k->k<0,A[iNeighbourEdges[1],:])[1]) # Since A[iNeighbourEdges[1],nextVertex[end]]>0, the next vertex downstream must have A[iNeighbourEdges[1],k]<0
#         else
#             # iNeighbourEdges[2] is downstream clockwise
#             push!(orderedEdges, iNeighbourEdges[2])
#             # Find the other vertex surrounding iNeighbourEdges[2] that isn't orderedVertices[end]
#             if A[iNeighbourEdges[2],nextVertex[end]] > 0 
#                 push!(orderedVertices, findall(k->k<0,A[iNeighbourEdges[2],:])[1])
#             else 
#                 push!(orderedVertices, findall(k->k>0,A[iNeighbourEdges[2],:])[1])
#             end
#         end 
#         nextVertex[1]=orderedVertices[end]
#     end 
    
#     # Convert to circular arrays 
#     orderedVertices = CircularArray(orderedVertices)
#     orderedEdges = CircularArray(orderedEdges)
    
#     return orderedVertices, orderedEdges

# end


# function orderAroundCell2(matrices, i)

#     @unpack A, B, C = matrices

#     # Find all edges and vertices for cell i
#     cellEdges, edgeOrientations = findnz(B[i,:])
#     cellVertices = findnz(C[i,:])[1]
    
#     orderedVertices = copy(cellVertices) # Ordered list of vertices around cell i in clockwise direction 
#     orderedEdges = copy(cellEdges)    # Ordered list of edges around cell i in clockwise direction
    
#     nextEdge = [1]
#     nextVertex = [1]
#     tst = [1]
#     for jj=1:length(cellEdges)
#         orderedEdges[jj] = cellEdges[nextEdge[1]]        
#         a, b = findnz(A[orderedEdges[jj],:])
#         b[1] == B[i,orderedEdges[jj]] ? tst[1] = a[1] : tst[1] = a[2]
#         nextVertex[1] = findfirst(x->x==tst[1],cellVertices)
#         orderedVertices[jj] = cellVertices[nextVertex[1]]
#         nextEdge[1] = findfirst(x->x==[x for x in cellEdges if x!=orderedEdges[jj] && A[x,orderedVertices[jj]]!=0][1],cellEdges)
#     end 
    
#     # Convert to circular arrays 
#     orderedVertices = CircularArray(orderedVertices)
#     orderedEdges = CircularArray(orderedEdges)
    
#     return orderedVertices, orderedEdges

# end

export orderAroundCell , orderAroundCell3 #, orderAroundCell3, orderAroundCell1

end