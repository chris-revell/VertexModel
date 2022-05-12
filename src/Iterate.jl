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
using FastBroadcast

# Local modules
# include("SpatialData.jl"); using .SpatialData
using SpatialData
# include("CalculateForce.jl"); using .CalculateForce
using CalculateForce
# include("T1Transitions.jl"); using .T1Transitions
using T1Transitions
# include("TopologyChange.jl"); using .TopologyChange
using TopologyChange
# include("Division.jl"); using .Division
using Division

function iterate!(iteration,params,matrices)

    @unpack R, tempR, ΔR, rkCoefficients, totalF = matrices
    @unpack dt,nCells = params

    totalF .= sum.(eachrow(matrices.F))
    @.. thread=false totalF .+= matrices.externalF
    @.. thread=false tempR .= R .+ totalF.*dt*rkCoefficients[1,iteration]
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

    totalF .= sum.(eachrow(matrices.F))
    @.. thread=false totalF .+= matrices.externalF
    @.. thread=false ΔR .+= totalF.*dt*rkCoefficients[2,iteration]/6.0

    return nothing
end

export iterate!

end
