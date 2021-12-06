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

function iterate!(iteration,rkCoefficients,R,tempR,ΔR,params,matrices)

    @unpack F = matrices
    @unpack dt = params

    tempR .= R .+ F.*dt*rkCoefficients[1,iteration]
    spatialData!(tempR,params,matrices)

    if iteration == 1
        if (t1Transitions!(R,params,matrices))==1
            topologyChange!(matrices)
            spatialData!(tempR,params,matrices)
        end
        fill!(ΔR,@SVector zeros(2))
    end
    calculateForce!(tempR,params,matrices)
    ΔR .+= F.*dt*rkCoefficients[2,iteration]/6.0

    return nothing
end



# function iterate2!(coefficient1,coefficient2,R,tempR,ΔR,params,matrices)
#
#     @unpack F = matrices
#     @unpack dt = params
#
#     tempR .= R .+ F.*dt/coefficient1
#     spatialData!(tempR,params,matrices)
#     calculateForce!(tempR,params,matrices)
#     ΔR .+= F.*dt/coefficient2
#
#     return nothing
# end

export iterate!#, iterate2!

end
