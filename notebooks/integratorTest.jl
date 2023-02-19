using DrWatson
using DifferentialEquations
using CairoMakie
using StaticArrays
using LinearAlgebra

dxLJ(r) = -(12*0.1^12)/r^13 + (6*0.1^6)/r^7

function model!(du,u,p,t)
    N = size(u)[2]
    fill!(du,0.0)
    for n=1:N
        for m in [a for a in 1:N if a!=n]
            x = norm(u[:,n].-u[:,m])            
            F = dxLJ(x).*(u[:,n].-u[:,m])./x            
            du[:,n] = du[:,n].-F
            du[:,m] = du[:,m].+F
        end
    end
    return du
end

u0 = rand(2,30)

prob = ODEProblem(model!,u0,(0.0,10.0),[])

# sol = solve(prob)

integrator = init(prob)

t=[0.0]

while integrator.t < 1.0
    
    step!(integrator)
    if 
    integrator.u = [integrator.u rand(2)]
end




fig = CairoMakie.Figure()
ax = Axis(fig[1,1],aspect=DataAspect())
xlims!(ax,(0,1))
ylims!(ax,(0,1))
    
mov = VideoStream(fig, framerate=5)

for u in sol.u
    empty!(ax)
    scatter!(ax,u[1,:],u[2,:],color=:black)
    recordframe!(mov)
end

save("test.mp4",mov)