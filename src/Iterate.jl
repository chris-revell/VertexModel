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

    @unpack R, tempR, ΔR, rkCoefficients = matrices
    @unpack dt = params

    tempR .= R .+ (sum(matrices.F,dims=2).+matrices.externalF).*dt*rkCoefficients[1,iteration]
    spatialData!(tempR,params,matrices)

    if iteration == 1

        if division!(params,matrices)>0
            topologyChange!(matrices)
            spatialData!(tempR,params,matrices)
        end
        if (t1Transitions!(tempR,params,matrices))>1
            topologyChange!(matrices)
            spatialData!(tempR,params,matrices)
        end

        fill!(ΔR,@SVector zeros(2))
    end

    calculateForce!(tempR,params,matrices)

    ΔR .+= (sum(matrices.F,dims=2).+matrices.externalF).*dt*rkCoefficients[2,iteration]/6.0

    return nothing
end

export iterate!

end
