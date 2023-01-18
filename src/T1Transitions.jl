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

function t1Transitions!(R,params,matrices,t)

    @unpack A,B,Ā,B̄,C,edgeLengths,edgeTangents,ϵ,boundaryVertices,edgeMidpoints,cellPositions = matrices
    @unpack nEdges,t1Threshold = params

    transitionCount = 0

    for j=1:nEdges
        if edgeLengths[j] < t1Threshold
            # Find vertices a and b at either end of the short edge j
            a = findall(j->j>0,A[j,:])[1]
            b = findall(j->j<0,A[j,:])[1]
            if boundaryVertices[a] == 0 && boundaryVertices[b] == 0 # Skip edges for which either vertex is at the boundary. Eventually we can probably figure out a better way of handling these edge cases                
                # Find cells around vertices a and b
                aCells = findall(i->i!=0,C[:,a])
                bCells = findall(i->i!=0,C[:,b])
                # Find cells P, Q, R, S surrounding vertices a and b
                Q = findall(i->i>0,B[:,j])[1] # Assume edge j has positive (clockwise) orientation with respect to cell Q
                S = findall(i->i<0,B[:,j])[1] # Assume edge j has negative (anti-clockwise) orientation with respect to cell S
                P = setdiff(aCells, [Q,S])[1] # Assume cell P shares vertex a, which has positive orientation with respect to edge j
                R = setdiff(bCells, [Q,S])[1] # Assume cell R shares vertex b, which has negative orientation
                # Find edges k and m to be removed from or added to vertices a and b 
                k = (findall(i->i!=0,B[S,:])∩findall(i->i!=0,B[P,:]))[1] # Edge k shared by cell P and cell S
                m = (findall(i->i!=0,B[Q,:])∩findall(i->i!=0,B[R,:]))[1] # Edge m shared by cell Q and cell R
                # Remove edge j from cells Q and S, assuming orientation from clockwise rotation of edge j
                B[Q,j] = 0; B[S,j] = 0
                # Add edge j to cells R and S, assuming orientation from clockwise rotation of edge j
                B[R,j] = 1; B[P,j] = -1
                # Add vertex b to edge k, setting orientation from matrix B
                A[k,b] = -B[P,k] #A[k,a]
                # Remove vertex a from edge k
                A[k,a] = 0
                # Add vertex a to edge m, setting orientation from matrix B
                A[m,a] = B[Q,m] #A[m,b]
                # Remove vertex b from edge m 
                A[m,b] = 0

                # # Change positions of vertices a and b. Ensure new edgeLengths[j] value is 10% longer than t1Threshold
                # R[a] = edgeMidpoints[j].+ϵ*edgeTangents[j].*1.1*t1Threshold/edgeLengths[j]
                # R[b] = edgeMidpoints[j].-ϵ*edgeTangents[j].*1.1*t1Threshold/edgeLengths[j]

                transitionCount += 1

                # Break loop when a T1 transition occurs, preventing more than 1 transition per time step. Eventually we can figure out a better way of handling multiple transitions per time step.
                break

            end
        end
    end

    return transitionCount

end

export t1Transitions!

end
