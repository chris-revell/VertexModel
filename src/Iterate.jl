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
using DrWatson
using FromFile

# Local modules
@from "SpatialData.jl" using SpatialData
@from "CalculateForce.jl" using CalculateForce
@from "T1Transitions.jl" using T1Transitions
@from "TopologyChange.jl" using TopologyChange
@from "Division.jl" using Division
@from "SenseCheck.jl" using SenseCheck

function iterate!(iteration,params,matrices,t)

    @unpack R, tempR, ΔR, rkCoefficients, totalF = matrices
    @unpack dt,nCells = params

    totalF .= sum.(eachrow(matrices.F))     # FastBroadcast doesn't work for this line; not sure why
    @.. thread=false totalF .+= matrices.externalF
    @.. thread=false tempR .= R .+ totalF.*dt*rkCoefficients[1,iteration]
    spatialData!(tempR,params,matrices)

    if iteration == 1

        if (t1Transitions!(tempR,params,matrices,t))>0
            senseCheck(matrices.A, matrices.B; marker="T1")
            topologyChange!(matrices)
            spatialData!(tempR,params,matrices)
        end
        if division!(params,matrices)>0
            senseCheck(matrices.A, matrices.B; marker="division")
            topologyChange!(matrices)
            spatialData!(tempR,params,matrices)
        end        
        fill!(ΔR,@SVector zeros(2))
    end

    calculateForce!(tempR,params,matrices)

    totalF .= sum.(eachrow(matrices.F))     # FastBroadcast doesn't work for this line; not sure why
    @.. thread=false totalF .+= matrices.externalF
    @.. thread=false ΔR .+= totalF.*dt*rkCoefficients[2,iteration]/6.0

    return nothing
end

export iterate!

end
