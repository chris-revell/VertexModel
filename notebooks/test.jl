# https://discourse.julialang.org/t/changing-size-of-ode-system-resize-works-but-deleteat-doesnt/45900/5

using DifferentialEquations
using Plots

const α = -0.3

function f(du,u,p,t)
    for i in 1:length(u)
        du[i] = α*u[i]
    end
end

function condition(u,t,integrator) # Event when event_f(u,t) == 0
    minimum(u)-0.01
end

function affect!(integrator)
#   u = integrator.u
  #resize!(integrator,length(u)-1) # This works
#   deleteat!(integrator, length(u)) # This errors
    integrator.u = [integrator.u; 6]
    nothing
end

callback = ContinuousCallback(condition,affect!)
u0 = [5, 3, 1]
tspan = (0.0,20.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob,Tsit5(),callback=callback)

Plots.plot(sol.t,map((x)->length(x),sol[:]),lw=3, ylabel="Number of Cells",xlabel="Time")