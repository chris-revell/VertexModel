#
#  T1Transitions.jl
#  VertexModelJL
#
#  Created by Christopher Revell on 10/03/2021.
#
#
# Function to recalculate derived matrices (Ā, Aᵀ etc.) for any change in vertex network topology (ie any change to A or B matrices).

module T1Transitions

# Julia packages
using LinearAlgebra

@inline @views function t1Transitions!(A,Ā,B,B̄,C,R,nEdges,edgeLengths,edgeTangents,t1Threshold,ϵ)

    transitionOccurred = 0

    for i=1:nEdges
        if edgeLengths[i] < t1Threshold

            # Find vertices a and b at either end of the short edge i
            a,b = findall(j->j!=0,Ā[i,:])

            # Find cells around vertices a and b
            aCells = findall(j->j!=0,C[:,a])
            bCells = findall(j->j!=0,C[:,b])

            # Find edges around vertices a and b
            aEdges = setdiff!(findall(j->j!=0,Ā[:,a]),[i])
            bEdges = setdiff!(findall(j->j!=0,Ā[:,b]),[i])

            # Find cells P,Q,R,S surrounding vertices a and b.
            # Cells Q and S are shared by both vertices a and b.
            # Cells R and P neighbour only cells a and b respectively.
            # Note that if vertex a or b is at the boundary, cell P or R will not exist.
            cellQ,cellS = intersect(aCells,bCells)
            cellR = setdiff(aCells,[cellQ,cellS])
            cellP = setdiff(bCells,[cellQ,cellS])

            # Set orientation for the whole system relative to cell Q
            originalEdgeOrientation = B[cellQ,i]

            # Remove edge i from cells Q and S
            B[cellQ,i] = 0
            B[cellS,i] = 0
            # Add edge i to cells R and S, being careful to avoid boundaries where there is no cell R or P
            length(cellR) > 0 ? B[cellR[1],i] = -originalEdgeOrientation : nothing
            length(cellP) > 0 ? B[cellP[1],i] = originalEdgeOrientation : nothing

            # Find edge z around cell Q and vertex a
            z = intersect(bEdges,findall(j->j!=0,B̄[cellQ,:]))[1]
            # Adjust A to reflect edge z now meeting vertex a, not b
            A[z,a] = A[z,b]
            A[z,b] = 0

            # Find edge e around cell S and vertex a
            y = intersect(aEdges,findall(j->j!=0,B̄[cellS,:]))[1]
            # Adjust A to reflect edge y now meeting vertex b, not a
            A[y,b] = A[y,a]
            A[y,a] = 0

            # Change positions of vertices a and b. Ensure new edgeLengths[i] value is 10% longer than t1Threshold
            R[a,:] .+= originalEdgeOrientation*0.55*t1Threshold.*ϵ*edgeTangents[i,:]./edgeLengths[i] .+ 0.5*edgeTangents[i,:]
            R[b,:] .-= originalEdgeOrientation*0.55*t1Threshold.*ϵ*edgeTangents[i,:]./edgeLengths[i] .+ 0.5*edgeTangents[i,:]

            transitionOccurred = 1
        end
    end

    return transitionOccurred

end

export t1Transitions!

end
