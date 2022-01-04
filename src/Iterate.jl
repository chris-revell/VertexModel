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
using StaticArrays

# Local modules
include("SpatialData.jl"); using .SpatialData
include("CalculateForce.jl"); using .CalculateForce
include("T1Transitions.jl"); using .T1Transitions
include("TopologyChange.jl"); using .TopologyChange
include("Division.jl"); using .Division

function iterate!(iteration,params,matrices)

    @unpack R, tempR, ΔR, F, rkCoefficients = matrices
    @unpack dt = params

    tempR .= R .+ F.*dt*rkCoefficients[1,iteration]
    spatialData!(tempR,params,matrices)

    if iteration == 1

        division!(params,matrices)
        if (t1Transitions!(tempR,params,matrices))==1
            topologyChange!(matrices)
            spatialData!(tempR,params,matrices)
        end

        fill!(ΔR,@SVector zeros(2))
    end

    #display(params.nCells)
    # display(matrices.R            )
    # display(matrices.tempR            )
    # display(matrices.ΔR               )
    # display(matrices.A                )
    # display(matrices.B                )
    # display(matrices.Aᵀ               )
    # display(matrices.Ā               )
    # display(matrices.Āᵀ              )
    # display(matrices.Bᵀ               )
    # display(matrices.B̄               )
    # display(matrices.B̄ᵀ              )
    # display(matrices.C                )
    # display(matrices.cellEdgeCount    )
    # display(matrices.boundaryVertices )
    # display(matrices.cellPositions    )
    # display(matrices.cellPerimeters   )
    # display(matrices.cellOrientedAreas)
    # display(matrices.cellAreas        )
    # display(matrices.cellTensions     )
    # display(matrices.cellPressures    )
    # display(matrices.cellAges         )
    # display(matrices.edgeLengths      )
    # display(matrices.edgeTangents     )
    # display(matrices.edgeMidpoints    )
    # display(matrices.vertexEdges      )
    # display(matrices.vertexCells      )
    # display(matrices.F                )

    calculateForce!(tempR,params,matrices)
    ΔR .+= F.*dt*rkCoefficients[2,iteration]/6.0

    return nothing
end

export iterate!

end
