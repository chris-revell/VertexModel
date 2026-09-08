#
#  Ablation.jl
#  VertexModel
#
#  Originally AblateCells.jl in Ablation.jl package
#

module Ablation

# Julia packages
using DiscreteCalculus
using SparseArrays

# Ablate all cells in list ablatedCellsList
function ablateCells(R, A, B, ablatedCellsList)
    newB = B[[ii for ii=1:size(B,1) if ii∉ablatedCellsList],[jj for jj=1:size(B,2)]]
    edgeCellPairs = [findall(x->x!=0,newB[:,j]) for j=1:size(newB,2)]
    orphanedEdges = findall(x->length(x)==0, edgeCellPairs)
    newB2 = newB[[ii for ii=1:size(newB,1)],[jj for jj=1:size(newB,2) if jj∉orphanedEdges]]
    newA = A[[jj for jj=1:size(A,1) if jj∉orphanedEdges],[kk for kk=1:size(A,2)]]
    edgeVertexPairs = [findall(x->x!=0,newA[:,k]) for k=1:size(newA,2)]
    orphanedVertices = findall(x->length(x)==0, edgeVertexPairs)
    newA2 = newA[[jj for jj=1:size(newA,1)],[kk for kk=1:size(newA,2) if kk∉orphanedVertices]]    
    # senseCheck(newA2, newB2)
    newR = [R[k] for k=1:length(R) if k∉orphanedVertices]
    return newR, newA2, newB2
end 

# Note, this function ablates a single edge j, unlike the ablateCells function
function ablateEdge(R, A, B, j)
    # Find cells adjacent to edge 
    j_is = findall(x->x!=0, tmpB[1][:,j])
    # Ablate adjacent cells 
    newR, newA, newB = ablateCells(R, A, B, j_is)
end
   
export ablateCells
export ablateEdge

end #end module 
