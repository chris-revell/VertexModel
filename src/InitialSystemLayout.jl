#
#  InitialSystemLayout.jl
#  VertexModel
#
#  Created by Christopher Revell on 06/09/2023.
#
#
# Function to create a hexagonal grid of cells. 
# Given number of rows nRows, central row has length nRows, each adjacent row has length nRows-1 etc. 
# Number of cells is then nRows*(nRows-1) - (floor(Int64, nRows/2)+1)*(floor(Int64, nRows/2)+2) + nRows

module InitialSystemLayout

# Julia packages
using LinearAlgebra
using SparseArrays
using StaticArrays
using DrWatson
using FromFile
using DelaunayTriangulation
using FromFile
using Random
using InvertedIndices
using NonlinearSolve
using Polynomials 

@from "SenseCheck.jl" using SenseCheck
@from "CubicSolutions.jl" using CubicSolutions

function initialEdgeLength(γ,L₀)
    a,b,c,d = 9/2, 0, (-√(3) + 12*γ), -2*γ*L₀
    p = Polynomial([d, c, b, a])
    roots_p = roots(p)
    # Only consider real roots
    tol = 1e-10
    real_roots = real.(roots_p[abs.(imag.(roots_p)) .< tol])
    if isempty(real_roots)
        error("No real roots from l cubic")
    else
        l = maximum(real_roots)
    end
    return l
end

function initialSystemLayout(γ,Λ_AA,Λ_BB,Λ_AB,Λ_AE,Λ_BE, nRows,spiky,initialSystem)

    # We start by assuming all cells are A-cells, grow the monolayer, and later assign B-cells
    L₀ = -Λ_AA/(2*γ)

    # plt = plot_parameter_space(100,Λ_AA,Λ_BB,γ)

    equilibriumEdgeLength = initialEdgeLength(γ, L₀)
    horizontalCellSpacing = 2.0*equilibriumEdgeLength*sin(π/3)
    verticalCellSpacing = 1.5*equilibriumEdgeLength

    t1Threshold = equilibriumEdgeLength*0.15

    l_AA = initialEdgeLength(γ, -Λ_AA/(2*γ))
    l_BB = initialEdgeLength(γ, -Λ_BB/(2*γ))
    l_AB = initialEdgeLength(γ, -Λ_AB/(2*γ))
    l_AE = initialEdgeLength(γ, -Λ_AE/(2*γ))
    l_BE = initialEdgeLength(γ, -Λ_BE/(2*γ))

    # Hexagonal packing constant: vertical offset between adjacent rows of unit-spaced hexagon centres
    Δy = sqrt(1 - 0.5^2)
 
    if initialSystem == "32-cell"
        # Diamond-shaped grid: central row of length nRows, tapering by one cell per row moving outwards
        cellPoints = [SVector(x, 0.0) for x = 1:nRows]
        for j = 1:(floor(Int64,nRows/2))
            for i = 1:nRows-j
                # Need to add a small amount of randomness to prevent errors in voronoi tessellation 
                push!(cellPoints, SVector(i + 0.5 * j + rand() * 0.001 - 0.0005, j * Δy + rand() * 0.001 - 0.0005))
                push!(cellPoints, SVector(i + 0.5 * j + rand() * 0.001 - 0.0005, -j * Δy + rand() * 0.001 - 0.0005))
            end
        end
    elseif initialSystem == "2-row"
        # Two full rows of hexagons, each of width nRows cells, offset by half a cell to give proper hex packing
        width = nRows
        cellPoints = SVector{2,Float64}[]
        for i = 1:width+2 # 2 extras as ghosts for tesselation cut-0ff
            # Bottom row
            push!(cellPoints, SVector(i + rand() * 0.001 - 0.0005, -Δy/2 + rand() * 0.001 - 0.0005))
        end
        for i = 1:width+2
            # Top row, offset by 0.5 in x for hexagonal packing
            push!(cellPoints, SVector(i + 0.5 + rand() * 0.001 - 0.0005, Δy/2 + rand() * 0.001 - 0.0005))
        end
        # Ghost rows above and below  for the tesselation cut-offs: 
        for i=1:width+2
            # Bottom ghost 
            push!(cellPoints, SVector(i -0.5+ rand() * 0.001 - 0.0005, -3*Δy/2 + rand() * 0.001 - 0.0005))
        end
        for i=1:width+2
            # Top ghost 
            push!(cellPoints, SVector(i +1+ rand() * 0.001 - 0.0005, 3*Δy/2 + rand() * 0.001 - 0.0005))
        end
    end


    # nRows = 9 # Must be an odd number
    
    xs = [x[1] for x in cellPoints]
    ys = [x[2] for x in cellPoints]

    ptsArray = zeros(Float64, (2, length(cellPoints)))
    for (i, point) in enumerate(cellPoints)
        ptsArray[1, i] = point[1]
        ptsArray[2, i] = point[2]
    end

    triangulation_unconstrained = triangulate(ptsArray)
    tessellation_constrained = voronoi(triangulation_unconstrained, clip=true)

    #Exclude points outside constraining boundary
    usableVertices = Int64[]
    for a in values(tessellation_constrained.polygons)
        push!(usableVertices, a...)
    end
    sort!(unique!(usableVertices))
    # outerVertices = setdiff(collect(1:num_polygon_vertices(tessellation_constrained)), usableVertices)

    # Map vertex indices in tessellation to vertex indices in incidence matrices (after excluding outer vertices)
    vertexIndexingMap = Dict(usableVertices .=> collect(1:length(usableVertices)))

    Rtmp = SVector.(tessellation_constrained.polygon_points[usableVertices])

    # Find pairs of vertices connected by edges in tessellation 
    # Use incidence matrix indexing for vertices, and exclude outer vertices 
    pairs = [(vertexIndexingMap[p[1]], vertexIndexingMap[p[2]]) for p in keys(tessellation_constrained.adjacent.adjacent) if p[1] ∈ usableVertices && p[2] ∈ usableVertices]
    # Ensure lowest index is first in tuple, and remove duplicates 
    orderedPairs = unique([(min(p...), max(p...)) for p in pairs])

    nVerts = length(Rtmp)
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

    
    if initialSystem == "2-row"
        # Deal with ghost cells
        ghostCells = vcat(
            1,
            width+2,
            width+3,
            2*width+4,
            (2*width + 5):(3*width + 6),
            (3*width + 7):(4*width + 8)
        )
        B = B[Not(ghostCells), :]

        ghostEdges = findall(j -> nnz(B[:,j])==0, axes(B,2))
        A = A[Not(ghostEdges), :]
        B = B[:, Not(ghostEdges)]

        ghostVertices = findall(k -> nnz(A[:,k])==0, axes(A,2))
        A = A[:, Not(ghostVertices)]
        Rtmp = Rtmp[Not(ghostVertices)]

        nCells = size(B, 1)
        nVerts = size(A, 2)
        nEdges = size(A, 1)

    elseif initialSystem == "32-cell"
        # Prune peripheral vertices with 2 edges that both belong to the same cell
        # Making the assumption that there will never be two such vertices adjacent to each other
        verticesToRemove = Int64[]
        edgesToRemove = Int64[]
        if !spiky
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
        end
        A = A[Not(edgesToRemove),Not(verticesToRemove)]
        B = B[:, Not(edgesToRemove)]
        Rtmp = Rtmp[Not(verticesToRemove)]
        senseCheck(A, B; marker="Error after removing peripheral vertices")
    end

    

    R = SVector{2, Float64}[]
    for r in Rtmp 
        push!(R, SVector(horizontalCellSpacing*(r[1] - (nRows-1)/2 - 1.0 ), horizontalCellSpacing*r[2]))
    end

    return A, B, R, t1Threshold, l_AA, l_BB, l_AB, l_AE, l_BE

end

export initialSystemLayout 

end
