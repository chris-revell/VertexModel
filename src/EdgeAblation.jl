#
#  EdgeAblation.jl
#  VertexModel
#
#
#
# Function that ablates a specified edge j, 
# removes adjacent cells, and removes trailing 
# vertices and edges left by ablation

module EdgeAblation

# Julia packages
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

function edgeAblation(j, params, matrices, integrator)
    
    @unpack A,
        B = matrices 
    @unpack nVerts, 
        nEdges, 
        nCells = params

    # Find cells adjacent to edge j 
    j_is = sort(findall(x->x!=0, B[:,j]))
    # Find trailing vertices adjacent to edge j left behind after edge ablation. These vertices will be pruned 
    j_ks = findall(x->x!=0, A[j,:])
    
    cellsToRemove = fill(true, nCells)
    edgesToRemove = fill(true, nEdges)
    vertsToRemove = fill(true, nVerts)
    uIndsToRemove = fill(true, nVerts*2) # Indices of vertex position x and y components in integrator.u
    cellsToRemove[j_is] .= false
    edgesToRemove[j] = false
    # Also find both other edges adjacent to pruned vertices, one of which will be removed for each vertex, whilst the other 
    # is extended to span the space left by the removed vertex
    for k in j_ks
        uIndsToRemove[2*k-1:2*k] .= false
        vertsToRemove[k] = false
        # Find edges j′ and j′′ surrounding k, excluding original edge j
        j′, j′′ = [jj for jj in findall(x->x!=0, A[:,k]) if jj!=j] 
        # Find second adjacent vertex k′ of edge j′
        k′ = [kk for kk in findall(x->x!=0, A[j′,:]) if kk!=k][1] 
        # Schedule edge j′ for removal
        edgesToRemove[j′] = false
        # Connect edge j′′ to k′, preserving direction of j′′ with respect to vertex to be removed, k
        A[j′′, k′] = A[j′′,k]
    end

    # Remove edges, vertices, and cells from matrices to create new A and B matrices 
    Btmp = B[cellsToRemove, edgesToRemove]
    Atmp = A[edgesToRemove, vertsToRemove]
    # Update stored incidence matrices in container object 
    matrices.A = Atmp
    matrices.B = Btmp

    # Resize most matrices in matrices container 
    resizeMatrices!(params, matrices, size(Atmp,2), size(Btmp,2), size(Btmp,1))
    # Some matrices need special treatment because their values cannot be inferred from A, B, and R, so we need to carefully delete specific values
    shrinkIndependentMatrices!(matrices, j_is, findall(x->false, edgesToRemove), j_ks)
    # Update stored number of cells, edges, and vertices
    params.nVerts = size(Atmp,2)
    params.nEdges = size(Atmp,1)
    params.nCells = size(Btmp,1)

    # Reduce size of domain in integrator 
    tmpu = integrator.u[uIndsToRemove]
    resize!(integrator, 2*size(Atmp,2))
    integrator.u .= tmpu
    u_modified!(integrator, true)
        
    return nothing 

end

export edgeAblation

end
