#
#  Callbacks.jl
#  VertexModel
#
#

module Callbacks

using OrdinaryDiffEq

using LinearAlgebra

# Steady state condition 
function conditionSteadyState(u, t, integrator)
    maximum(abs.(get_du(integrator))) < integrator.p[3] ? true : false
end

# Maximum time condition 
function conditiontMax(u, t, integrator)
    integrator.t <= integrator.p[1].tMax ? false : true
end

# function affectTerminate!(integrator)
#     # if conditionSteadyState() returns true, terminate integrator and pass successful return code
#     println("Terminate")
#     terminate!(integrator)
# end

export conditionSteadyState
export conditiontMax
# export affectTerminate!

end