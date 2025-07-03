#
#  Callbacks.jl
#  VertexModel
#
#  Created by Christopher Revell on 01/07/2025.
#
#

module Callbacks

using OrdinaryDiffEq
using DiffEqCallbacks
using LinearAlgebra

function conditionSteadyState(u, t, integrator)
    # @show maximum(norm.(get_du(integrator)))
    maximum(abs.(get_du(integrator))) < 100.0*integrator.opts.abstol ? true : false
    # if maximum(norm.(get_du(integrator))) < 100.0*integrator.opts.abstol 
    #     @show maximum(norm.(get_du(integrator)))
    #     return true
    # else
    #     return false
    # end
end
function conditiontMax(u, t, integrator)
    integrator.t <= integrator.p[1].tMax ? false : true
end
function affectTerminate!(integrator)
    # if conditionSteadyState() returns true, terminate integrator and pass successful return code
    println("Terminate")
    terminate!(integrator)
end
# function affectTerminateSS!(integrator)
#     # if conditionSteadyState() returns true, terminate integrator and pass successful return code
#     println("Terminate at steady state")
#     terminate!(integrator)
# end
# function affectTerminatetMax!(integrator)
#     # if conditionSteadyState() returns true, terminate integrator and pass successful return code
#     println("Terminate at tMax")
#     terminate!(integrator)
# end

# Create callback using user-defined functions above
# cbtMax = DiscreteCallback(conditiontMax, affectTerminate!)
# cbSS = DiscreteCallback(conditionSteadyState, affectTerminate!)

export conditionSteadyState
export conditiontMax
export affectTerminate!

end