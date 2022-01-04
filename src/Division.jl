#
#  Division.jl
#  VertexModel
#
#  Created by Christopher Revell on 13/12/2021.
#
#
# Function to divide cells that meet division criteria

module Division

# Julia packages
using LinearAlgebra
using StaticArrays
using SparseArrays
using UnPack

include("TopologyChange.jl"); using .TopologyChange

indexLoop(a,N) = (N+a-1)%(N)+1

function division!(params,matrices)

    @unpack nCells, nonDimCellCycleTime = params
    @unpack R, B, C, cellAges, edgeMidpoints, ϵ = matrices

    divisionCount = 0

    for i=1:nCells
        if cellAges[i] > nonDimCellCycleTime

            divisionCount += 1

            # Find all edges and vertices for cell i
            cellEdges = findall(j->j!=0,B[i,:])
            cellVertices = findall(j->j!=0,C[i,:])

            n = length(cellVertices)

            # Find long axis of cell by calculating the two furthest separated vertices
            maxDist = 0
            longPair = [0,0]
            for j=1:n-1
                for k=j+1:n
                    tmpDist = norm(R[cellVertices[j]].-R[cellVertices[k]])
                    if tmpDist>maxDist
                        maxDist = tmpDist
                        longPair .= [cellVertices[j],cellVertices[k]]
                    end
                end
            end
            # longAxis is the vector separating the two vertices with labels stored in longPair
            longAxis = R[longPair[1]].-R[longPair[2]]
            # centrePoint is the position halfway between the two vertices in longPair
            centrePoint = R[longPair[2]].+0.5.*longAxis

            # Find short axis perpendicular to long axis
            shortAxis = ϵ*longAxis
            shortAngle = (atan(shortAxis...)+2π)%2π

            # Calculate polar angles of vertices and edges around centrePoint
            vertexAngles = similar(cellVertices)
            edgeAngles = similar(cellEdges)
            for (k,v) in enumerate(cellVertices)
                vertexAngles[k] = atan(R[v]-centrePoint...)
            end
            for (k,e) in enumerate(cellEdges)
                edgeAngles[k]   = atan(edgeMidpoints[e]-centrePoint...)
            end

            # Adjust angles so that one vertex is at 0 and all other angles are relative
            minVertexAngle = minimum(vertexAngles)
            # Set vertex angles relative to minimum angle vertex
            vertexAngles .= (vertexAngles.-minVertexAngle).%2π
            # Use polar angles to sort cellVertices list
            cellVertices .= cellVertices[sortperm(vertexAngles)]
            sort!(vertexAngles)

            # Set edge angles relative to minimum angle vertex
            edgeAngles .= (edgeAngles.-minVertexAngle.+2π).%2π
            # Use polar angles to sort cellEdges list
            cellEdges .= cellEdges[sortperm(edgeAngles)]
            sort!(edgeAngles)

            # Find indices within cellEdges of edges intersected by short axis
            intersectedIndex = [n,n]
            for j=1:n-1
                if vertexAngles[j] <= shortAngle < vertexAngles[j+1]
                    intersectedIndex[1]= j
                elseif vertexAngles[j] <= (shortAngle+π)%2π < vertexAngles[j+1]
                    intersectedIndex[2]= j
                end
            end
            sort!(intersectedIndex)

            # Add new rows and columns to B matrix for new cell and edges
            B = [B; spzeros(1,nEdges)]
            B = [B spzeros(nCells+1,3)]

            # Add edges to new cell with existing orientations
            B[nCells+1,cellEdges[indexLoop(intersectedIndex[1]+1,n):intersectedIndex[2]-1]] .= B[i,cellEdges[indexLoop(intersectedIndex[1]+1,n):intersectedIndex[2]-1]]
            # Remove edges from existing cell that have been moved to new cell
            B[i,cellEdges[indexLoop(intersectedIndex[1]+1,n):intersectedIndex[2]-1]] .= 0
            # Add new short axis edge to existing cell
            B[i,nEdges+1] = 1
            # Add new short axis edge to new cell
            B[nCells+1,nEdges+1] = -1
            # Add new edges created by splitting intersected edges to new cell with the same orientation as the edge that was split
            B[nCells+1,nEdges+2] = B[i,cellEdges[intersectedIndex[1]]]
            B[nCells+1,nEdges+3] = B[i,cellEdges[intersectedIndex[2]]]

            # Find the two neighbouring cells that share the intersected edges
            neighbour1 = findall(j->j!=0,B[:,cellEdges[intersectedIndex[1]]])[1]
            neighbour2 = findall(j->j!=0,B[:,cellEdges[intersectedIndex[2]]])[1]

            # Add new edges to these neighbour cells
            B[neighbour1,nEdges+2] = B[neighbour1,cellEdges[intersectedIndex[1]]]
            B[neighbour2,nEdges+3] = B[neighbour2,cellEdges[intersectedIndex[2]]]

            # Add new rows and columns to A matrix for new vertices and edges
            A = [A; spzeros(3,nVerts)]
            A = [A spzeros(nEdges+3,2)]

            # Allocate new vertices to new and existing edges and existing vertices to new edges
            A[nEdges+2,cellVertices[indexLoop(intersectedIndex[1]+1,n)]] = A[cellEdges[intersectedIndex[1]],cellVertices[indexLoop(intersectedIndex[1]+1,n)]]
            A[nEdges+2,nVerts+1] = A[cellEdges[intersectedIndex[1]],cellVertices[intersectedIndex[1]]]
            A[cellEdges[intersectedIndex[1]],nVerts+1] = A[cellEdges[intersectedIndex[1]],cellVertices[indexLoop(intersectedIndex[1]+1,n)]]
            A[cellEdges[intersectedIndex[1]],cellVertices[indexLoop(intersectedIndex[1]+1,n)]] = 0

            A[nEdges+3,cellVertices[intersectedIndex[2]]] = A[cellEdges[intersectedIndex[2]],cellVertices[intersectedIndex[2]]]
            A[nEdges+3,nVerts+2] = A[cellEdges[intersectedIndex[2]],cellVertices[indexLoop(intersectedIndex[2]+1,n)]]
            A[cellEdges[intersectedIndex[2]],nVerts+2] = A[cellEdges[intersectedIndex[2]],cellVertices[intersectedIndex[2]]]
            A[cellEdges[intersectedIndex[2]],cellVertices[intersectedIndex[2]]] = 0

            A[nEdges+1,nVerts+1] = -1
            A[nEdges+1,nVerts+2] = 1

            # Add new vertex position
            newPos1 = (R[cellVertices[indexLoop(intersectedIndex[1]+1,n)]].+R[cellVertices[intersectedIndex[1]]])./2
            newPos2 = (R[cellVertices[indexLoop(intersectedIndex[2]+1,n)]].+R[cellVertices[intersectedIndex[2]]])./2
            R = vcat(R, [newPos1,newPos2])

        else
            #nothing
        end
    end


    if divisionCount>0

        params.nCells += 1*divisionCount
        params.nVerts += 2*divisionCount
        params.nEdges += 3*divisionCount

        # Add 1 component to vectors for new cell
        cellEdgeCount = vcat(cellEdgeCount,zeros(Int64,divisionCount))
        cellPositions = vcat(cellPositions,Array{SVector{2,Float64}}(undef,divisionCount))
        cellPerimeters = vcat(cellPerimeters,zeros(Float64,divisionCount))
        cellOrientedAreas = vcat(cellOrientedAreas,Array{SMatrix{2,2,Float64}}(undef,divisionCount))
        cellAreas = vcat(cellAreas,zeros(Float64,divisionCount))
        cellTensions = vcat(cellTensions,zeros(Float64,divisionCount))
        cellPressures = vcat(cellPressures,zeros(Float64,divisionCount))
        cellAges = vcat(cellAges,zeros(Float64,divisionCount))
        # Add 2 components to vectors for new vertices
        tempR = vcat(tempR,Array{SVector{2,Float64}}(undef,2*divisionCount))
        ΔR = vcat(ΔR,Array{SVector{2,Float64}}(undef,2*divisionCount))
        boundaryVertices = vcat(boundaryVertices,zeros(Int64,2*divisionCount))
        F = vcat(F,Array{SVector{2,Float64}}(undef,2*divisionCount))
        # Add 3 components to vectors for new edges
        edgeLengths = vcat(edgeLengths,zeros(Float64,3*divisionCount))
        edgeTangents = vcat(edgeTangents,Array{SVector{2,Float64}}(undef,3*divisionCount))
        edgeMidpoints = vcat(edgeMidpoints,Array{SVector{2,Float64}}(undef,3*divisionCount))


        Aᵀ = spzeros(Int64,nVerts,nEdges)
        Ā  = spzeros(Int64,nEdges,nVerts)
        Āᵀ = spzeros(Int64,nVerts,nEdges)
        Bᵀ = spzeros(Int64,nEdges,nCells)
        B̄  = spzeros(Int64,nCells,nEdges)
        B̄ᵀ = spzeros(Int64,nEdges,nCells)
        C  = spzeros(Int64,nCells,nVerts)
        topologyChange!(matrices)
    end

    return nothing

end

export division!

end
