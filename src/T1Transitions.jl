#
#  T1Transitions.jl
#  VertexModel
#
#  Created by Christopher Revell on 10/03/2021.
#
#
# 

module T1Transitions

# Julia packages
using LinearAlgebra
using UnPack
using StaticArrays

function t1Transitions!(integrator, params, matrices)

    @unpack A,
        B,
        C,
        edgeLengths,
        edgeTangents,
        timeSinceT1,
        boundaryEdges,
        boundaryCells,
        ϵ,
        firstT1onEdge,
        edgeLabels = matrices
    @unpack initialSystem,
        boundaryType,
        nEdges,
        t1Threshold,
        nonDimCycleTime,
        tMax,
        t1timeGap,
        l_AA,
        l_BB,
        l_AB,
        l_AE,
        l_BE,
        β = params

    # Reinterpret state vector as a vector of SVectors 
    R_u = reinterpret(SVector{2,Float64}, integrator.u)

    transitionCount = 0
    l_vec = [l_AA, l_AB, l_BB, l_AE, l_BE]
    t1ThresholdVec = zeros(nEdges)
    t1ThresholdVec .= l_vec[edgeLabels .+ 1]
    t1ThresholdVec .= t1ThresholdVec.*0.2

    for j=1:nEdges

        # if edgeLengths[j] < t1ThresholdVec[j] && timeSinceT1[j] > t1timeGap 
        #     println("t1 attempted but threshold time has not elapsed")
        # end
       
        
        if edgeLengths[j] < t1ThresholdVec[j] && (timeSinceT1[j] > t1timeGap || firstT1onEdge[j] == 0) 
        

            # Stopping the back and forth of t1s when the system equilibrates
            # if firstT1onEdge[j] >= 9
            #     println("T1 attempted and failed")
            #     break
            # end
            
            
            timeSinceT1[j] = 0

            # Find vertices a and b at either end of the short edge j
            a = findall(x -> x > 0, @view A[j, :])[1]
            b = findall(x -> x < 0, @view A[j, :])[1]
            # Find cells around vertices a and b
            aCells = findall(x -> x != 0, @view C[:, a])
            bCells = findall(x -> x != 0, @view C[:, b])
            # println("t1 transition triggerred at vertices ",a,",",b)
            

            if length(aCells) > 1 && length(bCells) > 1 # Exclude edges for which one vertex belongs to only one cell
                if boundaryType == "free" 
                    if boundaryEdges[j] == 0
                        # Find cells P, Q, R, S surrounding vertices a and b
                        Q = findall(x -> x > 0, @view B[:, j])[1] # Assume edge j has positive (clockwise) orientation with respect to cell Q
                        S = findall(x -> x < 0, @view B[:, j])[1] # Assume edge j has negative (anti-clockwise) orientation with respect to cell S                
                        aEdges = findall(x -> x != 0, @view A[:, a])                # Find all edges around vertex a
                        k = setdiff(aEdges, findall(x -> x != 0, @view B[Q, :]))[1]  # Find edge k around vertex a that is not shared by cell Q
                        bEdges = findall(x -> x != 0, @view A[:, b])                # Find all edges around vertex b
                        m = setdiff(bEdges, findall(x -> x != 0, @view B[S, :]))[1]  # Find edge m around vertex b that is not shared by cell S                    
                        # Assume cell P shares vertex a, which has positive orientation with respect to edge j
                        P = setdiff(aCells, [Q, S]) # NB This is an array that may have 1 element or be empty since cell P may not exist if vertex a is at the periphery, but the algorithm is generalised to accommodate this
                        # Assume cell R shares vertex b, which has negative orientation with respect to edge j
                        R = setdiff(bCells, [Q, S]) # NB This is an array that may have 1 element or be empty since cell R may not exist if vertex b is at the periphery, but the algorithm is generalised to accommodate this    
                        # Remove edge j from cells Q and S, assuming orientation from clockwise rotation of edge j
                        B[Q, j] = 0
                        B[S, j] = 0
                        # Add edge j to cells R and P, assuming orientation from clockwise rotation of edge j
                        B[R, j] .= 1
                        B[P, j] .= -1 # NB using . notation here and passing R and P as an array rather than a single value accommodates the possibility that R or P is an empty array
                        # Add vertex b to edge k, setting orientation from previous orientation of edge a
                        A[k, b] = A[k, a]
                        # Remove vertex a from edge k
                        A[k, a] = 0
                        # Add vertex a to edge m, setting orientation from previous orientation of edge b
                        A[m, a] = A[m, b]
                        # Remove vertex b from edge m 
                        A[m, b] = 0
                    else
                        
                        # Boundary edge 
                        # Find cells P, Q, R surrounding vertices a and b. There is no cell S.
                        Q = findall(x -> x != 0, @view B[:, j])[1]
                        # Assume cell P shares vertex a, which has positive orientation with respect to edge j
                        P = setdiff(aCells, [Q])[1]
                        # Assume cell R shares vertex b, which has negative orientation with respec to edge j
                        R = setdiff(bCells, [Q])[1]

                        aEdges = findall(x -> x != 0, @view A[:, a])                            # Find all edges around vertex a
                        k = setdiff(aEdges, findall(x -> x != 0, @view B[Q, :]))[1]              # Find edge k around vertex a that is not shared by cell Q
                        bEdges = findall(x -> x != 0, @view A[:, b])                            # Find all edges around vertex b
                        m = (findall(x -> x != 0, @view B[Q, :])∩findall(x -> x != 0, @view B[R, :]))[1]    # Find edge m around vertex b that is shared by Q and R

                        # Add edge j to cells R and P
                        B[R, j] = B[Q, j]
                        B[P, j] = -B[Q, j]
                        # Remove edge j from cell Q
                        B[Q, j] = 0
                        # Add vertex b to edge k, setting orientation from previous orientation of edge a
                        A[k, b] = A[k, a]
                        # Remove vertex a from edge k
                        A[k, a] = 0
                        # Add vertex a to edge m, setting orientation from previous orientation of edge b
                        A[m, a] = A[m, b]
                        # Remove vertex b from edge m 
                        A[m, b] = 0
                    end
                    R_u[b] = R_u[b] .+ 0.5.*edgeTangents[j] .+ 0.5.*ϵ*edgeTangents[j]
                    R_u[a] = R_u[a] .- 0.5.*edgeTangents[j] .- 0.5.*ϵ*edgeTangents[j]

                elseif boundaryType == "periodic"
                    # Find cells P, Q, R, S surrounding vertices a and b
                        Q = findall(x -> x > 0, @view B[:, j])[1] # Assume edge j has positive (clockwise) orientation with respect to cell Q
                        S = findall(x -> x < 0, @view B[:, j])[1] # Assume edge j has negative (anti-clockwise) orientation with respect to cell S                
                        aEdges = findall(x -> x != 0, @view A[:, a])                # Find all edges around vertex a
                        k = setdiff(aEdges, findall(x -> x != 0, @view B[Q, :]))[1]  # Find edge k around vertex a that is not shared by cell Q
                        bEdges = findall(x -> x != 0, @view A[:, b])                # Find all edges around vertex b
                        m = setdiff(bEdges, findall(x -> x != 0, @view B[S, :]))[1]  # Find edge m around vertex b that is not shared by cell S                    
                        # Assume cell P shares vertex a, which has positive orientation with respect to edge j
                        P = setdiff(aCells, [Q, S]) # NB This is an array that may have 1 element or be empty since cell P may not exist if vertex a is at the periphery, but the algorithm is generalised to accommodate this
                        # Assume cell R shares vertex b, which has negative orientation with respect to edge j
                        R = setdiff(bCells, [Q, S]) # NB This is an array that may have 1 element or be empty since cell R may not exist if vertex b is at the periphery, but the algorithm is generalised to accommodate this    
                        # Remove edge j from cells Q and S, assuming orientation from clockwise rotation of edge j
                        B[Q, j] = 0
                        B[S, j] = 0
                        # Add edge j to cells R and P, assuming orientation from clockwise rotation of edge j
                        B[R, j] .= 1
                        B[P, j] .= -1 # NB using . notation here and passing R and P as an array rather than a single value accommodates the possibility that R or P is an empty array
                        # Add vertex b to edge k, setting orientation from previous orientation of edge a
                        A[k, b] = A[k, a]
                        # Remove vertex a from edge k
                        A[k, a] = 0
                        # Add vertex a to edge m, setting orientation from previous orientation of edge b
                        A[m, a] = A[m, b]
                        # Remove vertex b from edge m 
                        A[m, b] = 0    

                    R_u[b] = R_u[b] .+ 0.5.*edgeTangents[j] .+ 0.5.*ϵ*edgeTangents[j]
                    R_u[a] = R_u[a] .- 0.5.*edgeTangents[j] .- 0.5.*ϵ*edgeTangents[j]

                end
                
                transitionCount += 1
                if firstT1onEdge[j] == 0
                    firstT1onEdge[j] = 1
                end
                
                # Break loop when a T1 transition occurs, preventing more than 1 transition per time step. Eventually we can figure out a better way of handling multiple transitions per time step.
                break
            elseif boundaryType == "free" && (length(aCells) == 1 || length(bCells) == 1)
                # The case where a vertex on j belongs only to one cell - i.e., a 'spike' at the periphery

                
                if length(bCells) == 1 # b is the spike, no cell R
                    println("first case t1 at spiky boundary edge j = ", j)
                    Q = findall(x -> x!=0, @view B[:,j])[1]
                    # Assume cell P shares vertex a
                    P = setdiff(aCells, [Q])[1]

                    aEdges = findall(x -> x!=0, @view A[:,a])
                    k = setdiff(aEdges, findall(x -> x!=0, @view B[Q,:]))[1] # Find edge k around vertex a that is not shared by cell Q
                    bEdges = findall(x -> x!=0, @view A[:,b])   
                    m = setdiff(bEdges,j)[1]# Find edge m around vertex b that IS shared by cell Q
                
                    # Add edge j to cell P, with the opposite orientation relative to cell Q
                    B[P, j] = -B[Q, j]
                    # Remove j from cell Q
                    B[Q,j] = 0
                    # Add vertex b to edge k, setting orientation from previous orientation of edge a
                    A[k, b] = A[k, a]
                    # Remove vertex a from edge k
                    A[k, a] = 0
                    # Add vertex a to edge m, setting orientation from previous orientation of edge b
                    A[m, a] = A[m, b]
                    # Remove vertex b from edge m 
                    A[m, b] = 0
                elseif length(aCells) == 1 # a is the spike 
                    println("second case t1 at spiky boundary edge j = ", j)
                    Q = findall(x -> x!=0, @view B[:,j])[1]
                    # Assume cell R shares vertex b
                    R = setdiff(bCells, [Q])[1]

                    bEdges = findall(x -> x!=0, @view A[:,b])   
                    m = (findall(x -> x != 0, @view B[Q, :])∩findall(x -> x != 0, @view B[R, :]))[1] # Find edge m around vertex b that IS shared by cell Q
                
                    # Add edge j to cell R, with the same orientation relative to cell Q
                    B[R, j] = B[Q, j]
                    # Remove j from cell Q
                    B[Q,j] = 0
                    # Add vertex a to edge m, setting orientation from previous orientation of edge b
                    A[m, a] = A[m, b]
                    # Remove vertex b from edge m 
                    A[m, b] = 0
                else
                    println("Spiky T1 failed - neither vertex is a spike.")
                    
                end 

                R_u[b] = R_u[b] .+ 0.5.*edgeTangents[j] .+ 0.5.*ϵ*edgeTangents[j]
                R_u[a] = R_u[a] .- 0.5.*edgeTangents[j] .- 0.5.*ϵ*edgeTangents[j]
                
                transitionCount += 1
                if firstT1onEdge[j] == 0
                    firstT1onEdge[j] = 1
                elseif β==0.0
                    firstT1onEdge[j] += 1
                    println(firstT1onEdge[j])
                end
                
                
                break 
            end 
        end
    end

    return transitionCount

end

export t1Transitions!

end
