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

        # if division!(params,matrices)>0
        #     topologyChange!(matrices)
        #     spatialData!(tempR,params,matrices)
        # end
        # if (t1Transitions!(tempR,params,matrices))>1
        #     topologyChange!(matrices)
        #     spatialData!(tempR,params,matrices)
        # end


        fill!(ΔR,@SVector zeros(2))
    end
    for i=1:params.nCells
        if isnan(R[i][1]) || isnan(R[i][2])
            throw()
        end
    end
    calculateForce!(tempR,params,matrices)
    ΔR .+= F.*dt*rkCoefficients[2,iteration]/6.0

    return nothing
end

export iterate!

end
