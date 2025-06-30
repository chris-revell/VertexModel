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
using OrdinaryDiffEq
using Distributions
using GeometryBasics
using Random
using DelimitedFiles

# Local modules
@from "OrderAroundCell.jl" using OrderAroundCell
@from "ResizeMatrices.jl" using ResizeMatrices
@from "Stretch.jl" using Stretch

function division!(integrator,params,matrices, folderName)

    @unpack nCells,
        nEdges,
        nVerts,
        nonDimCycleTime,
        realCycleTime,
        distLogNormal,
        L₀, γ, t1Threshold, tStretch,tMemChange, stretchType = params
    @unpack A, 
        B, C,
        cellTimeToDivide, 
        cellPositions, 
        edgeMidpoints, 
        cellEdgeCount, 
        cellPerimeters,
        cellVertexOrders, 
        cellEdgeOrders,
        cellShapeTensor, 
        boundaryEdges, 
        edgeTangents,
        ϵ, 
        μ, 
        Γ, R_membrane, Rt, R_final,
        cellLineage,
        cellGeneration,
        cellIndex = matrices

    # Reinterpret state vector as a vector of SVectors 
    R = reinterpret(SVector{2,Float64}, integrator.u)

    divisionCount = 0

    for i=1:nCells
        if cellTimeToDivide[i]<=0.0 && cellEdgeCount[i]>3 # Cell can only divide if it has more than 3 edges


            # Long and short axis from eigenvectors of shapetensor
            # Put some sort of tolerance that if eigenvalues are approx equal we randomly choose a division orientation, eg circ >0.95
            eigenVals, eigenVecs = eigen(matrices.cellShapeTensor[i]) # eigenvalues/vectors listed smallest to largest eigval.
            
            # circ = abs(eigenVals[1]/eigenVals[2]) # circularity
            # #for very circular cells randomly choose division axis
            # if circ > 0.9
            #     theta = rand()*π
            #     shortvec = [cos(theta), sin(theta)]
            # else
                if eigenVecs[:,1][2] < 0.0 #make it so vector is pointing in positive y direction (to fit with existing code in assigning new edges)
                    shortvec = -1.0*cellPerimeters[i].*eigenVecs[:,1] # Multiplication by cell perimeter ensures this axis is long enough to completely cross the cell; eigenvector has unit length otherwise
                else
                    shortvec = cellPerimeters[i].*eigenVecs[:,1] # Multiplication by cell perimeter ensures this axis is long enough to completely cross the cell; eigenvector has unit length otherwise
                end
            #end
            shortAxisLine = Line(Point{2,Float64}(matrices.cellPositions[i].+shortvec), Point{2,Float64}(matrices.cellPositions[i].-shortvec))

            # Test cell edges for an intersection
            poly = LineString(Point{2, Float64}.(R[cellVertexOrders[i][0:end]])) # Start and end with the same vertex by indexing circular array from 0 to end
            intersections = [intersects(line, shortAxisLine) for line in poly] #find which edges intersect and where
            intersectedIndices = findall(x->x!=0, first.(intersections))
            
            intersectedEdges = cellEdgeOrders[i][intersectedIndices]

            newCellVertices = cellVertexOrders[i][intersectedIndices[1]:intersectedIndices[2]-1] # Not including new vertices to be added later
            oldCellVertices = setdiff(cellVertexOrders[i], newCellVertices) # Not including new vertices to be added later
            newCellEdges = cellEdgeOrders[i][intersectedIndices[1]+1:intersectedIndices[2]-1] # IntersectedEdges allocated to old cell, not including new edge to be added later
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
                neighbourCell = setdiff(findall(x->x!=0, @view B[:,intersectedEdges[1]]),[i])[1]
                Btmp[neighbourCell,newEdges[2]] = -1
            end
            if boundaryEdges[intersectedEdges[2]] == 0
                neighbourCell = setdiff(findall(x->x!=0, @view B[:,intersectedEdges[2]]),[i])[1]
                Btmp[neighbourCell,newEdges[3]] = -1
            end

            # Ensure orientations with respect to other cells of all edges in new cell are correct
            for edge in [newCellEdges...]
                if boundaryEdges[edge] == 0
                    neighbourCell = setdiff(findall(x->x!=0, @view B[:,edge]),[i])[1]
                    Btmp[neighbourCell,edge] = -1
                end
            end

            # Add 3 new rows and 2 new columns to A matrix for new vertices and edges
            Atmp = spzeros(Int64,newEdges[end],newVertices[end])
            Atmp[1:nEdges,1:nVerts] .= A
            
            # First intersected old edge (which remains in the old cell) loses downstream vertex and gains the new vertex
            Atmp[intersectedEdges[1],cellVertexOrders[i][intersectedIndices[1]]] = 0
            Atmp[intersectedEdges[1],newVertices[1]] = A[intersectedEdges[1],cellVertexOrders[i][intersectedIndices[1]]]

            # Second intersected old edge (which remains in the old cell) loses upstream vertex and gains the new vertex             
            Atmp[intersectedEdges[2],cellVertexOrders[i][intersectedIndices[2]-1]] = 0
            Atmp[intersectedEdges[2],newVertices[2]] = A[intersectedEdges[2],cellVertexOrders[i][intersectedIndices[2]-1]]

            # First new edge gains both new vertices
            # Clockwise orientation of new edge with respect to new cell means the new edge leaves the second new vertex and enters the first new vertex
            Atmp[newEdges[1],newVertices[1]] = 1
            Atmp[newEdges[1],newVertices[2]] = -1

            # Second new edge gains first new vertex upstream (-1) and downstream vertex of old first intersected edge downstream (1)
            Atmp[newEdges[2],newVertices[1]] = -1
            Atmp[newEdges[2],cellVertexOrders[i][intersectedIndices[1]]] = 1
            
            # Third new edge gains second new vertex downstream (1) and upstream vertex of second old intersected edge upstream (-1)
            Atmp[newEdges[3],cellVertexOrders[i][intersectedIndices[2]-1]] = -1
            Atmp[newEdges[3],newVertices[2]] = 1

            # Check orientation of old vertices in new cell with respect to old edges allocated to new cell now that those old edges have been made clockwise with respect to new cell
            for (k,edge) in enumerate(newCellEdges)
                Atmp[edge, newCellVertices[k]] = -1
                Atmp[edge, newCellVertices[k+1]] = 1
            end

            # Add new vertex positions
            resize!(integrator,length(integrator.u)+4)
                        
            integrator.u[end-3:end-2] .= intersections[intersectedIndices[1]][2]
            integrator.u[end-1:end] .= intersections[intersectedIndices[2]][2]
            
            if stretchType !="none"
                params.tMemChange = integrator.t
                tStretch-integrator.t > 0 ? matrices.R_membrane .= matrices.Rt : matrices.R_membrane.=matrices.R_final
                
                # ccMem=(C*matrices.R_membrane)./cellEdgeCount 
                
                # shortAxisLineMem = Line(Point{2,Float64}(ccMem[i].+(2*shortvec)), Point{2,Float64}(ccMem[i].-(2*shortvec))) #division orientation should come from current cell shape
                    
                # polyMembrane = LineString(Point{2, Float64}.(matrices.R_membrane[cellVertexOrders[i][0:end]])) # Start and end with the same vertex by indexing circular array from 0 to end
                # intersectionsMembrane = [intersects(line, shortAxisLineMem) for line in polyMembrane] #might need to find shapetensor etc for R_membrane poly.
                # intersectedIndicesMem = findall(x->x!=0, first.(intersectionsMembrane))

                push!(matrices.R_membrane, intersections[intersectedIndices[1]][2])
                push!(matrices.R_membrane, intersections[intersectedIndices[2]][2])

                # #push!(matrices.R_membrane, intersectionsMembrane[intersectedIndicesMem[1]][2])
                # #push!(matrices.R_membrane, intersectionsMembrane[intersectedIndicesMem[2]][2])
                
            
                resize!(matrices.Rt,length(matrices.Rt)+2)
                resize!(matrices.R_final,length(matrices.R_final)+2)
                
            end


            matrices.A = Atmp
            matrices.B = Btmp
            resizeMatrices!(params, matrices, nVerts+2, nEdges+3, nCells+1)
            maxIndex=maximum(cellIndex)
            open(io -> writedlm(io, [integrator.t cellIndex[i] maxIndex+1 maxIndex+2 cellPositions[i][1] cellPositions[i][2] cellLineage[i] cellGeneration[i] eigenVecs[:,2][1] eigenVecs[:,2][2]], ','), datadir(folderName,"Division_$(@savename L₀ γ realCycleTime).csv"), "a") # write

            # Matrices not handled in resizeMatrices
            cellTimeToDivide[i] = rand(distLogNormal)*nonDimCycleTime
            cellGeneration[i]+=1
            cellIndex[i]=maxIndex+1
            push!(cellTimeToDivide,rand(distLogNormal)*nonDimCycleTime)
            push!(cellLineage,cellLineage[i])
            push!(cellGeneration,cellGeneration[i])
            push!(cellIndex, maxIndex+2)
            push!(matrices.μ, 1.0)
            push!(matrices.Γ, params.γ)

            #do we want a cell index so we can track daughter cells currently on split on cell retains parent's cell id

            divisionCount = 1
            break
        end 
    end
    return divisionCount

end

export division!

end
