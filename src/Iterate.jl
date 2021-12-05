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

function iterate1!(coefficient2,R,ΔR,params,matrices)

    @unpack F = matrices
    @unpack dt = params

    spatialData!(R,params,matrices)

    if (t1Transitions!(R,params,matrices))==1
        topologyChange!(matrices)
        spatialData!(R,params,matrices)
    end

    calculateForce!(R,params,matrices)
    ΔR .= F.*dt/coefficient2

    return nothing
end



function iterate2!(coefficient1,coefficient2,R,tempR,ΔR,params,matrices)

    @unpack F = matrices
    @unpack dt = params

    tempR .= R .+ F.*dt/coefficient1
    spatialData!(tempR,params,matrices)
    calculateForce!(tempR,params,matrices)
    ΔR .+= F.*dt/coefficient2

    return nothing
end

export iterate1!, iterate2!

end
