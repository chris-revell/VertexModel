#
#  Iterate.jl
#  VertexModel
#
#  Created by Christopher Revell on 03/12/2021.
#
#
# Function to perform one iteration of Runge-Kutta integration

module Iterate

# Julia packages
using LinearAlgebra
using UnPack

# Local modules
include("SpatialData.jl"); using .SpatialData
include("CalculateForce.jl"); using .CalculateForce
include("T1Transitions.jl"); using .T1Transitions
include("TopologyChange.jl"); using .TopologyChange

function iterate1!(coefficient2,R,ΔR,F,dt,nVerts,nCells,nEdges,γ,preferredPerimeter,preferredArea,A,B,Ā,B̄,C,cellEdgeCount,cellPositions,cellPerimeters,cellOrientedAreas,cellAreas,cellTensions,cellPressures,edgeLengths,edgeTangents,edgeMidpoints,boundaryVertices,ϵ,t1Threshold)

    # @unpack F = matrices
    # @unpack dt = params

    spatialData!(R,nCells,nEdges,γ,preferredPerimeter,preferredArea,A,B,Ā,B̄,C,cellEdgeCount,cellPositions,cellPerimeters,cellOrientedAreas,cellAreas,cellTensions,cellPressures,edgeLengths,edgeTangents,edgeMidpoints)

    if (t1Transitions!(R,nEdges,t1Threshold,A,B,Ā,B̄,C,edgeLengths,edgeTangents))==1
        topologyChange!(A,B,Aᵀ,Ā,Āᵀ,Bᵀ,B̄,B̄ᵀ,C,cellEdgeCount,boundaryVertices)
        spatialData!(R,nCells,nEdges,γ,preferredPerimeter,preferredArea,A,B,Ā,B̄,C,cellEdgeCount,cellPositions,cellPerimeters,cellOrientedAreas,cellAreas,cellTensions,cellPressures,edgeLengths,edgeTangents,edgeMidpoints)
    end

    calculateForce!(R,nVerts,nCells,nEdges,A,B,Ā,B̄,cellTensions,cellPressures,edgeLengths,edgeTangents,F,ϵ)
    ΔR .= F.*dt/coefficient2

    return nothing
end



function iterate2!(coefficient1,coefficient2,R,tempR,ΔR,F,dt,nVerts,nCells,nEdges,γ,preferredPerimeter,preferredArea,A,B,Ā,B̄,C,cellEdgeCount,cellPositions,cellPerimeters,cellOrientedAreas,cellAreas,cellTensions,cellPressures,edgeLengths,edgeTangents,edgeMidpoints,boundaryVertices,ϵ)

    # @unpack F = matrices
    # @unpack dt = params

    tempR .= R .+ F.*dt/coefficient1
    spatialData!(tempR,nCells,nEdges,γ,preferredPerimeter,preferredArea,A,B,Ā,B̄,C,cellEdgeCount,cellPositions,cellPerimeters,cellOrientedAreas,cellAreas,cellTensions,cellPressures,edgeLengths,edgeTangents,edgeMidpoints)
    calculateForce!(tempR,nVerts,nCells,nEdges,A,B,Ā,B̄,cellTensions,cellPressures,edgeLengths,edgeTangents,F,ϵ)
    ΔR .+= F.*dt/coefficient2

    return nothing
end

export iterate1!, iterate2!

end
