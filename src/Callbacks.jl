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

tss = TerminateSteadyState(abstol = 1e-8,
                    reltol = 1e-6,
                    test = allDerivPass; 
                    min_t = nothing,
                    wrap_test::Val = Val(true)
        )

function conditionSteadyState(u, t, integrator)
    maximum(norm.(get_du(integrator))) < 1e-3 ? true : false
end

function affectTerminate!(integrator)
    # if conditionSteadyState() returns true, terminate integrator and pass successful return code
    println("Terminate at steady state")
    terminate!(integrator, ReturnCode.Success)    
end
# Create callback using two user-defined functions above
cb = DiscreteCallback(conditionSteadyState, affectTerminate!)

