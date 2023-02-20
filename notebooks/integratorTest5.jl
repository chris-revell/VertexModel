using DrWatson
using DifferentialEquations
using CairoMakie
using StaticArrays
using LinearAlgebra
using GeometryBasics

N=10
tMax = 6.0

dxLJ(r) = -(12*0.1^12)/r^13 + (6*0.1^6)/r^7

function model!(du,u,p,t)
    fill!(du,@SVector zeros(2))
    for (n,i) in enumerate(u)
        for (m,j) in enumerate(u[n+1:end])
            F = dxLJ(norm(i-j)).*(i-j)/norm(i-j)
            du[n] -= F
            du[n+m] += F
        end
    end
    return du
end

function condition(u,t,integrator) # Event when event_f(u,t) == 0
    t%1.0 < integrator.dt ? (return 0) : (return 1)
end

# function affect!(integrator)
#     u = integrator.u
#     resize!(integrator,length(u)+1)
#     u[end] = @SVector rand(2)
#     nothing
# end

# callback = ContinuousCallback(condition,affect!)

u0 = [@SVector rand(2) for _=1:N]

tSpan = (0.,tMax)

prob = ODEProblem(model!, u0, tSpan)

integrator = init(prob, Tsit5())

fig = CairoMakie.Figure()
ax = Axis(fig[1,1],aspect=DataAspect())
xlims!(ax,(0,1))
ylims!(ax,(0,1))    
mov = VideoStream(fig, framerate=200)

while integrator.t < tMax
    empty!(ax)
    scatter!(ax,Point2.(first.(integrator.u),last.(integrator.u)),color=:red)
    recordframe!(mov)
    step!(integrator)
    if integrator.t%0.5 < integrator.dt
        # u = integrator.u
        resize!(integrator,length(integrator.u)+1)
        integrator.u[end] = @SVector rand(2)
        # push!(integrator.u,@SVector rand(2))
        u_modified!(integrator,true)
    end
end

save("test.mp4",mov)