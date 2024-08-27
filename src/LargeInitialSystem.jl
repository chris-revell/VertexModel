#
#  LargeInitialSystem.jl
#  VertexModel
#
#  Created by Christopher Revell on 06/09/2023.
#
#
# Function to create an initial system of 1 or 3 hexagonal cells.

module LargeInitialSystem

# Julia packages
using LinearAlgebra
using SparseArrays
using StaticArrays
using DrWatson
using FromFile
using DelaunayTriangulation
using FromFile
using Random

@from "SenseCheck.jl" using SenseCheck

function largeInitialSystem()

    cellPoints = [SVector(x, 0.0) for x=1:9]
    for j=1:4
        for i=1:9-j
            # Need to add a small amount of randomness to prevent errors in voronoi tessellation 
            push!(cellPoints, SVector(i + 0.5 * j + rand() * 0.001 - 0.0005, j * sqrt(1 - 0.5^2) + rand() * 0.001 - 0.0005))
            push!(cellPoints, SVector(i + 0.5 * j + rand() * 0.001 - 0.0005, -j * sqrt(1 - 0.5^2) + rand() * 0.001 - 0.0005))
        end
    end
    xs = [x[1] for x in cellPoints]
    ys = [x[2] for x in cellPoints]

    ptsArray = zeros(Float64, (2, length(cellPoints)))
    for (i, point) in enumerate(cellPoints)
        ptsArray[1, i] = point[1]
        ptsArray[2, i] = point[2]
    end

    triangulation_unconstrained = triangulate(ptsArray)
    tessellation_constrained = voronoi(triangulation_unconstrained, true)

    #Exclude points outside constraining boundary
    usableVertices = Int64[]
    for a in values(tessellation_constrained.polygons)
        push!(usableVertices, a...)
    end
    sort!(unique!(usableVertices))
    outerVertices = setdiff(collect(1:num_polygon_vertices(tessellation_constrained)), usableVertices)

    # Map vertex indices in tessellation to vertex indices in incidence matrices (after excluding outer vertices)
    vertexIndexingMap = Dict(usableVertices .=> collect(1:length(usableVertices)))

    R = SVector.(tessellation_constrained.polygon_points[usableVertices])

    # Find pairs of vertices connected by edges in tessellation 
    # Use incidence matrix indexing for vertices, and exclude outer vertices 
    pairs = [(vertexIndexingMap[p[1]], vertexIndexingMap[p[2]]) for p in keys(tessellation_constrained.adjacent.adjacent) if p[1] ∈ usableVertices && p[2] ∈ usableVertices]
    # Ensure lowest index is first in tuple, and remove duplicates 
    orderedPairs = unique([(min(p...), max(p...)) for p in pairs])

    nVerts = length(R)
    nEdges = length(orderedPairs)
    nCells = length(cellPoints)

    # Construct A matrix mapping tessellation edges to tessellation vertices 
    A = spzeros(Int64, nEdges, nVerts)
    for (edgeIndex, vertices) in enumerate(orderedPairs)
        A[edgeIndex, vertices[1]] = 1
        A[edgeIndex, vertices[2]] = -1
    end

    # NB get_polygon(tessellation_constrained,x) or tessellation_constrained.polygons[x] return indices of vertices around cell x ordered anti-clockwise, with first and last element the same

    # Construct B matrix mapping voronoi cell around each fibril to surrounding edges between vertices in tessellation
    # NB assume ϵᵢ is a clockwise rotation so cell orientation is into page. 
    B = spzeros(Int64, nCells, nEdges)
    for c = 1:nCells
        for i = 2:length(tessellation_constrained.polygons[c])
            vertexLeading = vertexIndexingMap[tessellation_constrained.polygons[c][i-1]]  # Leading with respect to *clockwise* direction around cell
            vertexTrailing = vertexIndexingMap[tessellation_constrained.polygons[c][i]]
            # Find index of edge connecting these vertices 
            edge = (findall(x -> x != 0, @view A[:, vertexLeading])∩findall(x -> x != 0, @view A[:, vertexTrailing]))[1]
            if A[edge, vertexLeading] > 0
                B[c, edge] = 1
            else
                B[c, edge] = -1
            end
        end
    end

#=

    # Prune peripheral vertices with 2 edges that both belong to the same cell
    # Making the assumption that there will never be two such vertices adjacent to each other
    verticesToRemove = Int64[]
    edgesToRemove = Int64[]
    for i = 1:nVerts
        edges = findall(x -> x != 0, @view A[:, i])
        cells1 = findall(x -> x != 0, @view B[:, edges[1]])
        cells2 = findall(x -> x != 0, @view B[:, edges[2]])
        if cells1 == cells2
            # If the lists of cells to which both edges of vertex i belong are identical, this implies that the edges are peripheral and only belong to one cell, so edge i should be removed.
            push!(verticesToRemove, i)
            push!(edgesToRemove, edges[1])
        end
    end
    for i in verticesToRemove
        edges = findall(x -> x != 0, @view A[:, i])
        otherVertexOnEdge1 = setdiff(findall(x -> x != 0, @view A[edges[1], :]), [i])[1]
        A[edges[2], otherVertexOnEdge1] = A[edges[2], i]
        A[edges[1], otherVertexOnEdge1] = 0
    end

    =#

    #remove external cells

    peripheralEdges=findall(x->x!=0,(ones(nCells)'*B)')
    borderEdges=findall(x->x!=0,vec(A*(ones(nVerts)'-abs.(ones(nCells)'*B)*abs.(A))'))

    verticesToRemove = [x[2] for x in findall(x->x!=0,abs.(ones(nCells)'*B)*abs.(A))]
    edgesToRemove = vcat(peripheralEdges, borderEdges)

    cellsToRemove=[]
    for j in peripheralEdges
        cells=findall(x->x!=0,B[:,j])
        push!(cellsToRemove,cells[1])
    end
    
    unique!(cellsToRemove)
    


 A = A[setdiff(1:size(A,1),edgesToRemove), setdiff(1:size(A,2),verticesToRemove)]
 B = B[setdiff(1:size(B,1),cellsToRemove), setdiff(1:size(B,2),edgesToRemove)]
 R = R[setdiff(1:size(R,1),verticesToRemove)]

    senseCheck(A, B; marker="Removing peropheral vertices")

    return A, B, R

end

export largeInitialSystem 

end
