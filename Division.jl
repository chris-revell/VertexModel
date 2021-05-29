#
#  Division.jl
#  VertexModelJL
#
#  Created by Christopher Revell on 18/02/2021.
#
#
# Function to update arrays for cell division for all cells that meet division conditions.

module Division

# Julia packages
using LinearAlgebra
using SparseArrays

# Local modules
using TopologyChange

function division!(A,B,C,R,cellAges,cellEdgeCount,cellPositions,edgeMidpoints,nonDimCycleTime,nCells,nEdges,nVerts)

    nCellsOld = nCells

    for i=1:nCellsOld
        # Cell divides if its age exceeds the cell cycle time
        if cellAges[i] > nonDimCycleTime

            cellVertices = findall(j->j!=0,C[i,:])
            cellEdges    = findall(j->j!=0,B[i,:])

            # Add new column to matrix B to accommodate the new edge
            newBColumn = spzeros(nCells,1)
            B = [B newBColumn]

            # Separate routines for cells with odd and even vertex counts
            if iseven(cellEdgeCount[i])

                # Find longest axis of cell by calculating distances between all vertex pairs
                dxMax = 0.0
                maxPair = [0,0]
                for (a,aVert) in enumerate(cellVertices)
                    for (b,bVert) in enumerate(cellVertices)
                        dx = norm(R[aVert,:].-R[bVert,:])
                        if dx > dxMax
                            dxMax = dx
                            maxPair .= [aVert,bVert]
                        end
                    end
                end

                # Find polar angles of the vertices between which a new edge is added
                aPolarAngle = atan(R[maxPair[1],1] .- cellPositions[i,1], R[maxPair[1],2] .- cellPositions[i,2])
                bPolarAngle = atan(R[maxPair[2],1] .- cellPositions[i,1], R[maxPair[2],2] .- cellPositions[i,2])

                # Add a new edge between the two vertices with the greatest distance between them (i and j) to matrix A
                newARow = spzeros(1,nVerts)
                if aPolarAngle>bPolarAngle
                    newARow[1,maxPair[1]] = 1
                    newARow[1,maxPair[2]] = -1
                else
                    newARow[1,maxPair[1]] = -1
                    newARow[1,maxPair[2]] = 1
                end
                A = [A;newARow]
                nEdges = nEdges+1

                # Find edges around cell and polar angles of corresponding edgeMidpoints
                newBRow = spzeros(1,nEdges)
                for (j,k) in enumerate(cellEdges)
                    edgePolarAngle = atan(edgeMidpoints[k,1] .- cellPositions[i,1], edgeMidpoints[k,2] .- cellPositions[i,2])
                    if min(aPolarAngle,bPolarAngle) < edgePolarAngle < max(aPolarAngle,bPolarAngle)
                        newBRow[1,k] = B[i,k]
                        B[i,k] = 0
                    end
                end
                nCells = nCells+1
                B = [B;newBRow]
                B[i,nEdges] = 1
                B[nCells,nEdges] = -1

                cellAges[i] = 0.0

            else


                # Find longest axis of cell by calculating distances between all vertex pairs and opposing edge midpoints
                edgePolarAngles = Array{Float64}[]
                for (j,k) in enumerate(cellEdges)
                    push!(edgePolarAngles,atan(edgeMidpoints[k,1] .- cellPositions[i,1], edgeMidpoints[k,2] .- cellPositions[i,2]))


                dxMax = 0.0
                maxPair = [0,0]
                for (a,aVert) in enumerate(cellVertices)
                    for (b,bVert) in enumerate(cellVertices)
                        dx = norm(R[aVert,:].-R[bVert,:])
                        if dx > dxMax
                            dxMax = dx
                            maxPair .= [aVert,bVert]
                        end
                    end
                end




            end
        end
    end

    return A, B, nCells, nEdges

end

export division!

end
