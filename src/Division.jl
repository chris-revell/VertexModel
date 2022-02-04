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

#include("TopologyChange.jl"); using .TopologyChange

indexLoop(a,N) = (N+a-1)%(N)+1

function division!(params,matrices)

    @unpack nCells, nEdges, nVerts, nonDimCycleTime = params
    @unpack R, A, B, C, cellAges, edgeMidpoints, cellEdgeCount, cellPositions, cellPerimeters, cellOrientedAreas, cellAreas, cellTensions, cellPressures, tempR, ΔR, boundaryVertices, F, edgeLengths, edgeTangents, ϵ = matrices

    divisionCount = 0

    Atmp = copy(A)
    Btmp = copy(B)
    newRs = Array{SVector{2,Float64}}(undef,0)
    nCellsLocal = nCells
    nEdgesLocal = nEdges
    nVertsLocal = nVerts

    for i=1:nCells
        if cellAges[i] > nonDimCycleTime

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
            vertexAngles = zeros(size(cellVertices))
            edgeAngles = zeros(size(cellEdges))
            for (k,v) in enumerate(cellVertices)
                vertexAngles[k] = atan((R[v].-centrePoint)...)
            end
            for (k,e) in enumerate(cellEdges)
                edgeAngles[k]   = atan((edgeMidpoints[e].-centrePoint)...)
            end

            # Adjust angles so that one vertex is at 0 and all other angles are relative
            minVertexAngle = minimum(vertexAngles)
            # Set vertex angles relative to minimum angle vertex
            vertexAngles .= (vertexAngles.-minVertexAngle).%2π
            shortAngle = (shortAngle-minVertexAngle+2π)%2π

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
            Btmp = [Btmp; spzeros(1,nEdgesLocal)]
            Btmp = [Btmp spzeros(nCellsLocal+1,3)]

            # Add edges to new cell with existing orientations
            Btmp[nCellsLocal+1,cellEdges[indexLoop(intersectedIndex[1]+1,n):intersectedIndex[2]-1]] .= Btmp[i,cellEdges[indexLoop(intersectedIndex[1]+1,n):intersectedIndex[2]-1]]
            # Remove edges from existing cell that have been moved to new cell
            Btmp[i,cellEdges[indexLoop(intersectedIndex[1]+1,n):intersectedIndex[2]-1]] .= 0
            # Add new short axis edge to existing cell
            Btmp[i,nEdgesLocal+1] = 1
            # Add new short axis edge to new cell
            Btmp[nCellsLocal+1,nEdgesLocal+1] = -1
            # Add new edges created by splitting intersected edges to new cell with the same orientation as the edge that was split
            Btmp[nCellsLocal+1,nEdgesLocal+2] = Btmp[i,cellEdges[intersectedIndex[1]]]
            Btmp[nCellsLocal+1,nEdgesLocal+3] = Btmp[i,cellEdges[intersectedIndex[2]]]

            # Find the neighbouring cells that share the intersected edges
            # Add new edges to these neighbour cells
            neighbours1 = symdiff(findall(j->j!=0,B[:,cellEdges[intersectedIndex[1]]]),[i])
            if size(neighbours1)[1] > 0
                Btmp[neighbours1[1],nEdgesLocal+2] = Btmp[neighbours1[1],cellEdges[intersectedIndex[1]]]
            end
            neighbours2 = symdiff(findall(j->j!=0,B[:,cellEdges[intersectedIndex[2]]]),[i])
            if size(neighbours2)[1] > 0
                Btmp[neighbours2[1],nEdgesLocal+3] = Btmp[neighbours2[1],cellEdges[intersectedIndex[2]]]
            end

            # Add new rows and columns to A matrix for new vertices and edges
            Atmp = [Atmp; spzeros(3,nVertsLocal)]
            Atmp = [Atmp spzeros(nEdgesLocal+3,2)]

            # Allocate new vertices to new and existing edges and existing vertices to new edges
            Atmp[nEdgesLocal+2,cellVertices[indexLoop(intersectedIndex[1]+1,n)]] = Atmp[cellEdges[intersectedIndex[1]],cellVertices[indexLoop(intersectedIndex[1]+1,n)]]
            Atmp[nEdgesLocal+2,nVertsLocal+1] = Atmp[cellEdges[intersectedIndex[1]],cellVertices[intersectedIndex[1]]]
            Atmp[cellEdges[intersectedIndex[1]],nVertsLocal+1] = Atmp[cellEdges[intersectedIndex[1]],cellVertices[indexLoop(intersectedIndex[1]+1,n)]]
            Atmp[cellEdges[intersectedIndex[1]],cellVertices[indexLoop(intersectedIndex[1]+1,n)]] = 0

            Atmp[nEdgesLocal+3,cellVertices[intersectedIndex[2]]] = Atmp[cellEdges[intersectedIndex[2]],cellVertices[intersectedIndex[2]]]
            Atmp[nEdgesLocal+3,nVertsLocal+2] = Atmp[cellEdges[intersectedIndex[2]],cellVertices[indexLoop(intersectedIndex[2]+1,n)]]
            Atmp[cellEdges[intersectedIndex[2]],nVertsLocal+2] = Atmp[cellEdges[intersectedIndex[2]],cellVertices[intersectedIndex[2]]]
            Atmp[cellEdges[intersectedIndex[2]],cellVertices[intersectedIndex[2]]] = 0

            Atmp[nEdgesLocal+1,nVertsLocal+1] = -1
            Atmp[nEdgesLocal+1,nVertsLocal+2] = 1

            # Add new vertex position
            newPos1 = (R[cellVertices[indexLoop(intersectedIndex[1]+1,n)]].+R[cellVertices[intersectedIndex[1]]])./2
            newPos2 = (R[cellVertices[indexLoop(intersectedIndex[2]+1,n)]].+R[cellVertices[intersectedIndex[2]]])./2
            # Add new vertex positions at axis intersection with existing edges using line intersection https://en.wikipedia.org/wiki/Line–line_intersection
            # a = edgeTangents[intersectedIndex[1]][2]/edgeTangents[intersectedIndex[1]][1]
            # c = R[intersectedIndex[1]][2] - R[intersectedIndex[1]][1]*a
            # b = shortAxis[2]/shortAxis[1]
            # d = centrePoint[2] - centrePoint[1]*b
            # newPos1 = SVector{2}([(d-c)/(a-b), a*(d-c)/(a-b)+c])
            # a = edgeTangents[intersectedIndex[2]][2]/edgeTangents[intersectedIndex[2]][1]
            # c = R[intersectedIndex[2]][2] - R[intersectedIndex[2]][1]*a
            # b = shortAxis[2]/shortAxis[1]
            # d = centrePoint[2] - centrePoint[1]*b
            # newPos2 = SVector{2}([(d-c)/(a-b), a*(d-c)/(a-b)+c])
            append!(newRs,[newPos1,newPos2])

            cellAges[i] = 0

            nCellsLocal += 1
            nVertsLocal += 2
            nEdgesLocal += 3

            break

        else
            #nothing
        end
    end

    if divisionCount>0

        params.nCells = nCellsLocal
        params.nVerts = nVertsLocal
        params.nEdges = nEdgesLocal

        # Add 1 component to vectors for new cell
        append!(cellEdgeCount,zeros(Int64,divisionCount))
        append!(cellPositions,Array{SVector{2,Float64}}(undef,divisionCount))
        append!(cellPerimeters,zeros(Float64,divisionCount))
        append!(cellOrientedAreas,Array{SMatrix{2,2,Float64}}(undef,divisionCount))
        append!(cellAreas,zeros(Float64,divisionCount))
        append!(cellTensions,zeros(Float64,divisionCount))
        append!(cellPressures,zeros(Float64,divisionCount))
        append!(cellAges,zeros(Float64,divisionCount))

        # Add 2 components to vectors for new vertices
        append!(R,newRs)
        append!(tempR,newRs)
        append!(ΔR,Array{SVector{2,Float64}}(undef,2*divisionCount))
        append!(boundaryVertices,zeros(Int64,2*divisionCount))
        # append!(F,Matrix{SVector{2,Float64}}(undef,2*divisionCount))
        # Add 3 components to vectors for new edges
        append!(edgeLengths,zeros(Float64,3*divisionCount))
        append!(edgeTangents,Array{SVector{2,Float64}}(undef,3*divisionCount))
        append!(edgeMidpoints,Array{SVector{2,Float64}}(undef,3*divisionCount))

        matrices.A = Atmp
        matrices.B = Btmp
        matrices.Aᵀ = spzeros(Int64,nVertsLocal,nEdgesLocal)
        matrices.Ā  = spzeros(Int64,nEdgesLocal,nVertsLocal)
        matrices.Āᵀ = spzeros(Int64,nVertsLocal,nEdgesLocal)
        matrices.Bᵀ = spzeros(Int64,nEdgesLocal,nCellsLocal)
        matrices.B̄  = spzeros(Int64,nCellsLocal,nEdgesLocal)
        matrices.B̄ᵀ = spzeros(Int64,nEdgesLocal,nCellsLocal)
        matrices.C  = spzeros(Int64,nCellsLocal,nVertsLocal)        
        matrices.F  = Matrix{SVector{2,Float64}}(undef,nVertsLocal,nCellsLocal)
        #fill!(F,@SVector zeros(2))

    end

    return divisionCount

end

export division!

end
