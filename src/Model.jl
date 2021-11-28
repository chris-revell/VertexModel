#
#  Model.jl
#  VertexModelJL
#
#  Created by Christopher Revell on 31/01/2021.
#
#

module Model

# Julia packages

# Local modules
include("TopologyChange.jl"); using .TopologyChange
include("CreateRunDirectory.jl"); using .CreateRunDirectory
include("SpatialData.jl"); using .SpatialData
include("CalculateForce.jl"); using .CalculateForce
include("T1Transitions.jl"); using .T1Transitions

function model(du, u, p, t)

    A,Ā,Aᵀ,Āᵀ,B,B̄,Bᵀ,B̄ᵀ,C,R,F,externalF,nCells,nVerts,nEdges,cellPositions,cellEdgeCount,cellAreas,cellOrientedAreas,cellPerimeters,cellTensions,cellPressures,edgeLengths,edgeMidpoints,edgeTangents,boundaryVertices,vertexEdges,gamma,preferredPerimeter,preferredArea,t1Threshold,ϵ,pressureExternal = p

    spatialData!(A,Ā,B,B̄,C,R,nCells,nEdges,cellPositions,cellEdgeCount,cellAreas,cellOrientedAreas,cellPerimeters,cellTensions,cellPressures,edgeLengths,edgeMidpoints,edgeTangents,gamma,preferredPerimeter,preferredArea)
    transitionOccurred = t1Transitions!(A,Ā,B,B̄,C,R,nEdges,edgeLengths,edgeTangents,t1Threshold,ϵ)
    if transitionOccurred==1
        topologyChange!(A,Ā,Aᵀ,Āᵀ,B,B̄,Bᵀ,B̄ᵀ,C,R,cellEdgeCount,cellPositions,boundaryVertices,vertexEdges,edgeTangents,nVerts,nCells)
        spatialData!(A,Ā,B,B̄,C,R,nCells,nEdges,cellPositions,cellEdgeCount,cellAreas,cellOrientedAreas,cellPerimeters,cellTensions,cellPressures,edgeLengths,edgeMidpoints,edgeTangents,gamma,preferredPerimeter,preferredArea)
        transitionOccurred = 0
    end
    calculateForce!(F,externalF,A,Ā,B,B̄,cellPressures,cellTensions,edgeTangents,edgeLengths,nVerts,nCells,nEdges,ϵ,pressureExternal,boundaryVertices)

    ΔR .= (F .+ externalF).*dt/6.0



    end






    outputToggle==1 ? run(`convert -delay 0 -loop 0 output/$folderName/plot"*".png output/$folderName/animated.gif`) : nothing

end

export vertexModel

end
