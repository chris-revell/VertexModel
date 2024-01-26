#
#  InitialHexagons.jl
#  VertexModel
#
#  Created by Natasha Cowley on 05/01/2024.
#
#
# Function to create an initial large random system

module RandomInitialSystem

# Julia packages
using LinearAlgebra
using SparseArrays
using StaticArrays
using DrWatson
using FromFile
using DelaunayTriangulation
using StaticArrays
using FromFile
using Random

@from "SenseCheck.jl" using SenseCheck

function randomInitialSystem()

    rng=MersenneTwister(1240)

    cellPoints = [SVector(x, 0.0) for x=1:15]
    for j=1:7
        for i=1:15-j
            # Need to add a small amount of randomness to prevent errors in voronoi tessellation 
            push!(cellPoints,SVector(i+0.5*j+(rand(rng)-0.5)*0.5, j*sqrt(1-0.5^2)+(rand(rng)-0.5)*0.5))
            push!(cellPoints,SVector(i+0.5*j+(rand(rng)-0.5)*0.5, -j*sqrt(1-0.5^2)+(rand(rng)-0.5)*0.5))
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
        push!(usableVertices,a...)
    end
    sort!(unique!(usableVertices))
    outerVertices = setdiff(collect(1:num_polygon_vertices(tessellation_constrained)),usableVertices)

    # Map vertex indices in tessellation to vertex indices in incidence matrices (after excluding outer vertices)
    vertexIndexingMap = Dict(usableVertices.=>collect(1:length(usableVertices)))

    R = SVector.(tessellation_constrained.polygon_points[usableVertices])

    # Find pairs of vertices connected by edges in tessellation 
    # Use incidence matrix indexing for vertices, and exclude outer vertices 
    pairs = [(vertexIndexingMap[p[1]],vertexIndexingMap[p[2]]) for p in keys(tessellation_constrained.adjacent.adjacent) if p[1]∈usableVertices && p[2]∈usableVertices]
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
            edge = (findall(x -> x != 0, A[:, vertexLeading])∩findall(x -> x != 0, A[:, vertexTrailing]))[1]
            if A[edge, vertexLeading] > 0
                B[c, edge] = 1
            else
                B[c, edge] = -1
            end
        end
    end

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

#Random T1s

    nCells=length(B[:,1])
    nEdges=length(B[1,:])
    nVerts=length(R)    

    C=0.5.*abs.(B)*abs.(A)
tangents=A*R
e_len=norm.(tangents)

for _=1:1
    peripheralEdges=findall(x->x!=0,(ones(nCells)'*B)')
    borderEdges=findall(x->x!=0,vec(A*(ones(nVerts)'-abs.(ones(nCells)'*B)*abs.(A))'))
    interiorEdges=setdiff(1:nEdges, vcat(peripheralEdges, borderEdges))
    shortEdges=intersect(findall(x->x<0.35,e_len), interiorEdges)
    
    randomEdges=randsubseq(rng, interiorEdges, 0.05)    
    for j in randomEdges
            
            a,b=findall(x->x!=0,A[j, :])


            aCells = findall(i->i!=0,C[:,a])            
            bCells = findall(i->i!=0,C[:,b])            
            if length(aCells)>1 && length(bCells) > 1 # Exclude edges for which one vertex belongs to only one cell
                    # Find cells P, Q, R, S surrounding vertices a and b
                    Q = findall(i->i>0,B[:,j])[1] # Assume edge j has positive (clockwise) orientation with respect to cell Q
                    S = findall(i->i<0,B[:,j])[1] # Assume edge j has negative (anti-clockwise) orientation with respect to cell S                
                    aEdges = findall(x->x!=0,A[:,a])                # Find all edges around vertex a
                    k = setdiff(aEdges,findall(x->x!=0,B[Q,:]))[1]  # Find edge k around vertex a that is not shared by cell Q
                    bEdges = findall(x->x!=0,A[:,b])                # Find all edges around vertex b
                    m = setdiff(bEdges,findall(x->x!=0,B[S,:]))[1]  # Find edge m around vertex b that is not shared by cell S                    
                    # Assume cell P shares vertex a, which has positive orientation with respect to edge j
                    P = setdiff(aCells, [Q,S]) # NB This is an array that may have 1 element or be empty since cell P may not exist if vertex a is at the periphery, but the algorithm is generalised to accommodate this
                    # Assume cell R shares vertex b, which has negative orientation with respect to edge j
                    K = setdiff(bCells, [Q,S]) # NB This is an array that may have 1 element or be empty since cell R may not exist if vertex b is at the periphery, but the algorithm is generalised to accommodate this    
                    # Remove edge j from cells Q and S, assuming orientation from clockwise rotation of edge j
                    B[Q,j] = 0; B[S,j] = 0
                    # Add edge j to cells R and P, assuming orientation from clockwise rotation of edge j
                    B[K,j] .= 1; B[P,j] .= -1 # NB using . notation here and passing R and P as an array rather than a single value accommodates the possibility that R or P is an empty array
                    # Add vertex b to edge k, setting orientation from previous orientation of edge a
                    A[k,b] = A[k,a]
                    # Remove vertex a from edge k
                    A[k,a] = 0
                    # Add vertex a to edge m, setting orientation from previous orientation of edge b
                    A[m,a] = A[m,b]
                    # Remove vertex b from edge m 
                    A[m,b] = 0
                
            end
            C=0.5.*abs.(B)*abs.(A)
            tangents=A*R
            e_len=norm.(tangents)
    end
end

senseCheck(A, B; marker="Random T1s")

    
    return A, B, R

end

export randomInitialSystem 

end
