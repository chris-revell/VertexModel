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
using DrWatson
using FromFile
using UnPack
using LinearAlgebra
using StaticArrays
using SparseArrays
using CircularArrays
using DifferentialEquations
using Distributions
using Random

# Local modules
@from "OrderAroundCell.jl" using OrderAroundCell
@from "ResizeMatrices.jl" using ResizeMatrices

function division!(integrator,params,matrices)

    @unpack nCells,
        nEdges,
        nVerts,
        nonDimCycleTime,
        distLogNormal,
        γ = params
    @unpack A, 
        B, 
        cellTimeToDivide, 
        cellPositions, 
        edgeMidpoints, 
        cellEdgeCount, 
        cellVertexOrders, 
        cellEdgeOrders, 
        boundaryEdges, 
        ϵ, 
        μ, 
        Γ = matrices

    divisionCount = 0

    newRs = Array{SVector{2,Float64}}(undef,0) # Positions of new vertices created by division

    for i=1:nCells
        # if μ[i]<1.5 && cellTimeToDivide[i]<=0.0 && cellEdgeCount[i]>3 # Cell can only divide if it has more than 3 edges
        if cellTimeToDivide[i]<=0.0 && cellEdgeCount[i]>3 # Cell can only divide if it has more than 3 edges

            # Find long axis of cell by calculating the two furthest separated vertices
            distances = zeros(Float64,1)
            longPair = [1,3]
            # TODO this block sometimes fails to find longPair indices, causing an index of 0 to be passed to longAxis = integrator.u[longPair[1]].-integrator.u[longPair[2]] 
            for j=1:cellEdgeCount[i]-1
                for k=j+1:cellEdgeCount[i]
                    tmpDist = norm(integrator.u[cellVertexOrders[i][j]].-integrator.u[cellVertexOrders[i][k]])
                    if tmpDist>maximum(distances)
                        longPair .= [cellVertexOrders[i][j],cellVertexOrders[i][k]]
                    end
                    push!(distances,tmpDist)
                end
            end
            # longAxis is the vector separating the two vertices with labels stored in longPair
            longAxis = integrator.u[longPair[1]].-integrator.u[longPair[2]]
            # Find short axis perpendicular to long axis
            shortAxis = ϵ*longAxis
            shortAngle1 = (atan(shortAxis...)+2π)%2π
            # Find the indices within cellEdgeOrders[i] of the edges intersected by the short axis
            intersectedIndex = [0,0]
            minAngle = atan((edgeMidpoints[cellEdgeOrders[i][1]].-cellPositions[i])...)
            for ind = 1:cellEdgeCount[i]
                if shortAngle1 <= atan((integrator.u[cellVertexOrders[i][ind]].-cellPositions[i])...)-minAngle
                    intersectedIndex[1] = ind
                    break
                end
            end
            intersectedIndex[2] = intersectedIndex[1]+cellEdgeCount[i]÷2
            intersectedEdges = cellEdgeOrders[i][intersectedIndex]

            newCellVertices = cellVertexOrders[i][intersectedIndex[1]:intersectedIndex[2]-1] # Not including new vertices to be added later
            oldCellVertices = setdiff(cellVertexOrders[i], newCellVertices) # Not including new vertices to be added later
            newCellEdges = cellEdgeOrders[i][intersectedIndex[1]+1:intersectedIndex[2]-1] # IntersectedEdges allocated to old cell, not including new edge to be added later
            oldCellEdges = setdiff(cellEdgeOrders[i], newCellEdges) # Not including new edge to be added later

            # Labels of new edges, cell, and vertices 
            newEdges = nEdges.+collect(1:3)
            newCell = nCells +1
            newVertices = nVerts.+collect(1:2)

            # Add 1 new row and 3 new columns to B matrix for new cell and 3 new edges
            Btmp = spzeros(Int64,newCell,newEdges[end])
            Btmp[1:nCells ,1:nEdges] .= B            

            # Add edges to new cell with clockwise orientations
            Btmp[newCell,newCellEdges] .= 1
            # Remove edges from existing cell that have been moved to new cell
            Btmp[i,newCellEdges] .= 0
            # Add new edge dividing cells to existing cell with anticlockwise orientation
            Btmp[i,newEdges[1]] = -1
            # Add all new edges to new cell with clockwise orientation
            Btmp[newCell,newEdges] .= 1
            
            # Find the neighbouring cells that share the intersected edges
            # Add new edges to these neighbour cells
            # These edges are clockwise in the new cell, so must be anticlockwise in these neighbour cells                      
            if boundaryEdges[intersectedEdges[1]] == 0
                neighbourCell = setdiff(findall(j->j!=0, @view B[:,intersectedEdges[1]]),[i])[1]
                Btmp[neighbourCell,newEdges[2]] = -1
            end
            if boundaryEdges[intersectedEdges[2]] == 0
                neighbourCell = setdiff(findall(j->j!=0, @view B[:,intersectedEdges[2]]),[i])[1]
                Btmp[neighbourCell,newEdges[3]] = -1
            end

            # Ensure orientations with respect to other cells of all edges in new cell are correct
            for edge in [newCellEdges...]
                if boundaryEdges[edge] == 0
                    neighbourCell = setdiff(findall(j->j!=0, @view B[:,edge]),[i])[1]
                    Btmp[neighbourCell,edge] = -1
                end
            end

            # Add 3 new rows and 2 new columns to A matrix for new vertices and edges
            Atmp = spzeros(Int64,newEdges[end],newVertices[end])
            Atmp[1:nEdges,1:nVerts] .= A
            
            # First intersected old edge (which remains in the old cell) loses downstream vertex and gains the new vertex
            Atmp[intersectedEdges[1],cellVertexOrders[i][intersectedIndex[1]]] = 0
            Atmp[intersectedEdges[1],newVertices[1]] = A[intersectedEdges[1],cellVertexOrders[i][intersectedIndex[1]]]

            # Second intersected old edge (which remains in the old cell) loses upstream vertex and gains the new vertex             
            Atmp[intersectedEdges[2],cellVertexOrders[i][intersectedIndex[2]-1]] = 0
            Atmp[intersectedEdges[2],newVertices[2]] = A[intersectedEdges[2],cellVertexOrders[i][intersectedIndex[2]-1]]

            # First new edge gains both new vertices
            # Clockwise orientation of new edge with respect to new cell means the new edge leaves the second new vertex and enters the first new vertex
            Atmp[newEdges[1],newVertices[1]] = 1
            Atmp[newEdges[1],newVertices[2]] = -1

            # Second new edge gains first new vertex upstream (-1) and downstream vertex of old first intersected edge downstream (1)
            Atmp[newEdges[2],newVertices[1]] = -1
            Atmp[newEdges[2],cellVertexOrders[i][intersectedIndex[1]]] = 1
            
            # Third new edge gains second new vertex downstream (1) and upstream vertex of second old intersected edge upstream (-1)
            Atmp[newEdges[3],cellVertexOrders[i][intersectedIndex[2]-1]] = -1
            Atmp[newEdges[3],newVertices[2]] = 1

            # Check orientation of old vertices in new cell with respect to old edges allocated to new cell now that those old edges have been made clockwise with respect to new cell
            for (k,edge) in enumerate(newCellEdges)
                Atmp[edge, newCellVertices[k]] = -1
                Atmp[edge, newCellVertices[k+1]] = 1
            end

            # Add new vertex positions
            resize!(integrator,length(integrator.u)+2)
            integrator.u[end-1:end] .= [edgeMidpoints[intersectedEdges[1]],edgeMidpoints[intersectedEdges[2]]]

            matrices.A = Atmp
            matrices.B = Btmp
            resizeMatrices!(params, matrices, nVerts+2, nEdges+3, nCells +1)
            # Matrices not handled in resizeMatrices
            cellTimeToDivide[i] = rand(distLogNormal)*nonDimCycleTime
            push!(cellTimeToDivide,rand(distLogNormal)*nonDimCycleTime)
            push!(matrices.μ, 1.0)
            push!(matrices.Γ, params.γ)
            push!(matrices.A₀s, params.A₀)

            divisionCount = 1
            
            break

        end 
    end

    return divisionCount

end

export division!

end
