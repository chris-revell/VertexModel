using DrWatson
using DifferentialEquations
using CairoMakie
using StaticArrays
using LinearAlgebra
using GeometryBasics

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

u0 = [@SVector rand(2) for _=1:N]

tSpan = (0.,10.)

prob = ODEProblem(model!, u0, tSpan)

sol = solve(prob, Tsit5())

fig = CairoMakie.Figure()
ax = Axis(fig[1,1],aspect=DataAspect())
xlims!(ax,(0,1))
ylims!(ax,(0,1))
    
mov = VideoStream(fig, framerate=100)

for u in sol.u
    empty!(ax)
    scatter!(ax,Point2.(first.(u),last.(u)),color=:black)
    recordframe!(mov)
end

save("test.mp4",mov)