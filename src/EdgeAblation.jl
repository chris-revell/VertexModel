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
using DifferentialEquations

@from "SenseCheck.jl" using SenseCheck
@from "ResizeMatrices.jl" using ResizeMatrices

function edgeAblation(j, params, matrices, integrator)
    
    @unpack A, B = matrices 
    @unpack nVerts, nEdges, nCells = params

    # Find cells adjacent to edge j 
    cells = sort(findall(x->x!=0, B[:,j]))

    # cells[1] will be a hole; cells[2] will be deleted
    # Set stiffness of cells[1] to 0
    matrices.μ[cells[1]] = 0.0
    matrices.Γ[cells[1]] = 0.0

    # Find edges in cell 2 that will be transferred to cell 1
    # No need to consider vertices transfered to cell 1 because these come free with adjustment to A matrix
    cell2EdgesTransferredToCell1 = sort([jj for jj in findall(x->x!=0, B[cells[2],:]) if jj!=j])

    # Find trailing vertices adjacent to edge j left behind after edge ablation. These vertices will be pruned 
    verticesToRemove = findall(x->x!=0,A[j,:])
    # Also find both other edges adjacent to pruned vertices, one of which will be removed for each vertex
    edgesToRemove = [j]
    for i in verticesToRemove
        # Edges around i excluding j
        edges = [jj for jj in findall(x->x!=0,A[:,i]) if jj!=j] 
        # Select first edge from edges around i, then find its other adjacent vertex, otherVertexOnEdge1
        otherVertexOnEdge1 = [kk for kk in findall(x->x!=0, A[edges[1],:]) if kk!=i][1] 
        # Attach the second edge from edges around i to the other vertex on first edge; preserve orientation of second edge        
        A[edges[2], otherVertexOnEdge1] = A[edges[2],i]
        # Unlink second edge from vertex otherVertexOnEdge1
        A[edges[2], i] = 0
        # Add edge to list of edges to remove 
        push!(edgesToRemove,edges[1])
        # Unlink the same edge from its other adjacent cell, unless the ablated edge is adjacent to a peripheral edge 
        if matrices.boundaryVertices[otherVertexOnEdge1] != 0
            otherCellToRemoveEdgeFrom = [ii for ii in findall(x->x!=0,B[:,edges[1]]) if ii∉cells][1]
            B[otherCellToRemoveEdgeFrom,edges[1]] = 0
        end
    end
    
    sort!(edgesToRemove)
    filter!(x->x∉edgesToRemove,cell2EdgesTransferredToCell1)

    # Adjust the indices of edges transferred to cell 1 to accommodate reduced size of new matrices. 
    # The same edge might have a lower index in the new A or B matrices than in the original matrices due to other deleted edges with lower indices.
    cell2EdgesTransferredToCell1Adjusted = Int64[]
    for x=1:length(cell2EdgesTransferredToCell1)
        if cell2EdgesTransferredToCell1[x] < edgesToRemove[1]
            push!(cell2EdgesTransferredToCell1Adjusted,cell2EdgesTransferredToCell1[x])
        elseif cell2EdgesTransferredToCell1[x] < edgesToRemove[2]
            push!(cell2EdgesTransferredToCell1Adjusted,cell2EdgesTransferredToCell1[x]-1)
        elseif cell2EdgesTransferredToCell1[x] < edgesToRemove[3]
            push!(cell2EdgesTransferredToCell1Adjusted,cell2EdgesTransferredToCell1[x]-2)
        else 
            push!(cell2EdgesTransferredToCell1Adjusted,cell2EdgesTransferredToCell1[x]-3)
        end
    end

    # Remove edges, vertices, and one cell from matrices to create new A and B matrices 
    newB = B[[ii for ii=1:nCells if ii!=cells[2]],[jj for jj=1:nEdges if jj∉edgesToRemove]]
    newB[cells[1],cell2EdgesTransferredToCell1Adjusted] .= B[cells[2],cell2EdgesTransferredToCell1]
    newA = A[[jj for jj=1:nEdges if jj∉edgesToRemove],[ii for ii=1:nVerts if ii∉verticesToRemove]]
    
    # Update stored incidence matrices in container object 
    matrices.A = newA
    matrices.B = newB
    
    resizeMatrices!(params, matrices, nVerts-2, nEdges-3, nCells-1)
    # Some matrices need special treatment because their values cannot be inferred from A, B, and R, so we need to delete specific values
    deleteat!(matrices.cellAges, cells[2])
    deleteat!(matrices.μ, cells[2])
    deleteat!(matrices.Γ, cells[2])

    # Reduce size of domain in integrator 
    deleteat!(integrator,verticesToRemove)
    u_modified!(integrator,true)
        
    return nothing 

end

export edgeAblation

end
