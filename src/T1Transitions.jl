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

function t1Transitions!(R,params,matrices)

    @unpack A,B,Ā,B̄,C,edgeLengths,edgeTangents,timeSinceT1,ϵ,boundaryVertices,boundaryEdges,edgeMidpoints,cellPositions = matrices
    @unpack nEdges,t1Threshold = params

    transitionCount = 0

    for j=1:nEdges
        if edgeLengths[j] < t1Threshold && timeSinceT1[j] > 1.0

            timeSinceT1[j] = 0

            # Find vertices a and b at either end of the short edge j
            a = findall(j->j>0,A[j,:])[1]
            b = findall(j->j<0,A[j,:])[1]
            # Find cells around vertices a and b
            aCells = findall(i->i!=0,C[:,a])            
            bCells = findall(i->i!=0,C[:,b])            
            if length(aCells)>1 && length(bCells) > 1 # Exclude edges for which one vertex belongs to only one cell
                if boundaryEdges[j]==0 
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
                    R = setdiff(bCells, [Q,S]) # NB This is an array that may have 1 element or be empty since cell R may not exist if vertex b is at the periphery, but the algorithm is generalised to accommodate this    
                    # Remove edge j from cells Q and S, assuming orientation from clockwise rotation of edge j
                    B[Q,j] = 0; B[S,j] = 0
                    # Add edge j to cells R and P, assuming orientation from clockwise rotation of edge j
                    B[R,j] .= 1; B[P,j] .= -1 # NB using . notation here and passing R and P as an array rather than a single value accommodates the possibility that R or P is an empty array
                    # Add vertex b to edge k, setting orientation from previous orientation of edge a
                    A[k,b] = A[k,a]
                    # Remove vertex a from edge k
                    A[k,a] = 0
                    # Add vertex a to edge m, setting orientation from previous orientation of edge b
                    A[m,a] = A[m,b]
                    # Remove vertex b from edge m 
                    A[m,b] = 0
                else
                    # Boundary edge 
                    # Find cells P, Q, R surrounding vertices a and b. There is no cell S.
                    Q = findall(i->i!=0,B[:,j])[1]
                    # Assume cell P shares vertex a, which has positive orientation with respect to edge j
                    P = setdiff(aCells, [Q])[1]   
                    # Assume cell R shares vertex b, which has negative orientation with respec to edge j
                    R = setdiff(bCells, [Q])[1]                                    
            
                    aEdges = findall(x->x!=0,A[:,a])                            # Find all edges around vertex a
                    k = setdiff(aEdges,findall(x->x!=0,B[Q,:]))[1]              # Find edge k around vertex a that is not shared by cell Q
                    bEdges = findall(x->x!=0,A[:,b])                            # Find all edges around vertex b
                    m = (findall(x->x!=0,B[Q,:])∩findall(x->x!=0,B[R,:]))[1]    # Find edge m around vertex b that is shared by Q and R
                    
                    # Add edge j to cells R and P
                    B[R,j] = B[Q,j]; B[P,j] = -B[Q,j]
                    # Remove edge j from cell Q
                    B[Q,j] = 0                
                    # Add vertex b to edge k, setting orientation from previous orientation of edge a
                    A[k,b] = A[k,a]
                    # Remove vertex a from edge k
                    A[k,a] = 0
                    # Add vertex a to edge m, setting orientation from previous orientation of edge b
                    A[m,a] = A[m,b]
                    # Remove vertex b from edge m 
                    A[m,b] = 0
                end
                transitionCount += 1
                # Break loop when a T1 transition occurs, preventing more than 1 transition per time step. Eventually we can figure out a better way of handling multiple transitions per time step.
                break               
            end 
        end
    end
    @show transitionCount
    return transitionCount

end

export t1Transitions!

end
