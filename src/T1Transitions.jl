#
#  T1Transitions.jl
#  VertexModel
#
#  Created by Christopher Revell on 10/03/2021.
#
#
# Function to recalculate derived matrices (Ā, Aᵀ etc.) for any change in vertex network topology (ie any change to A or B matrices).

module T1Transitions

# Julia packages
using LinearAlgebra
using UnPack

@views function t1Transitions!(R,params,matrices)

    @unpack A,B,Ā,B̄,C,edgeLengths,edgeTangents,ϵ,boundaryVertices = matrices
    @unpack nEdges,t1Threshold = params

    transitionCount = 0

    for j=1:nEdges
        if edgeLengths[j] < t1Threshold

            # Find vertices a and b at either end of the short edge j
            a = findall(j->j>0,A[j,:])
            b = findall(j->j<0,A[j,:])

            if boundaryVertices[a] != 0 || boundaryVertices[b] != 0
                # Skip edges for which either vertex is at the boundary
                # Eventually we can probably figure out a better way of handling these edge cases
            else
                println("T1 occurs")
                # Find cells around vertices a and b
                aCells = findall(i->i!=0,C[:,a])
                bCells = findall(i->i!=0,C[:,b])

                # Find edges around vertices a and b, not including j
                aEdges = setdiff!(findall(j->j!=0,A[:,a]),[j])
                aEdgesAngles = [
                    (atan((edgeMidpoints[aEdges[1]].-R[a])...)-atan((edgeMidpoints[j].-R[a])...)+2π)%2π,
                    (atan((edgeMidpoints[aEdges[2]].-R[a])...)-atan((edgeMidpoints[j].-R[a])...)+2π)%2π
                ]
                k,l = aEdges[sortperm(aEdgesAngles)]

                bEdges = setdiff!(findall(j->j!=0,A[:,b]),[j])
                bEdgesAngles = [
                    (atan((edgeMidpoints[bEdges[1]].-R[b])...)-atan((edgeMidpoints[j].-R[b])...)+2π)%2π,
                    (atan((edgeMidpoints[bEdges[2]].-R[b])...)-atan((edgeMidpoints[j].-R[b])...)+2π)%2π
                ]
                m,n = bEdges[sortperm(bEdgesAngles)]

                cellP = setdiff(aCells,bCells)[1]
                bCellAngles = [
                    (atan((cellPositions[bCells[1]].-R[b])...)-atan((edgeMidpoints[j].-R[b])...)+2π)%2π,
                    (atan((cellPositions[bCells[2]].-R[b])...)-atan((edgeMidpoints[j].-R[b])...)+2π)%2π,
                    (atan((cellPositions[bCells[3]].-R[b])...)-atan((edgeMidpoints[j].-R[b])...)+2π)%2π
                ]
                cellQ,cellR,cellS = bCells[sortperm(bCellAngles)]

                # Remove edge j from cells Q and S
                B[cellQ,j] = 0
                B[cellS,j] = 0
                # Add edge j to cells R and S
                B[cellR,j] = 1
                B[cellP,j] = -1

                A[k,b] = A[k,a]
                A[k,a] = 0

                A[m,a] = A[m,b]
                A[m,b] = 0

                # Change positions of vertices a and b. Ensure new edgeLengths[j] value is 10% longer than t1Threshold
                R[a] = edgeMidpoints[j].+ϵ*edgeTangents[j]./1.9
                R[b] = edgeMidpoints[j].-ϵ*edgeTangents[j]./1.9

                transitionCount += 1

                # Break loop when a T1 transition occurs
                # preventing more than 1 transition per time step.
                # Eventually we can figure out a better way of
                # handling multiple transitions per time step
                break

            end
        end
    end

    return transitionCount

end

export t1Transitions!

end
