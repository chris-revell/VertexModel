using CairoMakie
using StaticArrays
using Random
using Distributions
using LinearAlgebra
using GeometryBasics

nFibres = 100
lFibre = 300.0

xRange = (0.0,1000.0)
yRange = (0.0,1000.0)

fig = Figure()
ax = Axis(fig[1,1])

#%%

allFibres = Line[]
for i=1:nFibres
    pos = rand(Uniform(xRange...),2)
    orientation = normalize!(rand(Uniform(-1.0,1.0),2))
    line = Line(Point2(pos.-lFibre.*orientation./2), Point2(pos.+lFibre.*orientation./2))
    push!(allFibres,line)
end

#%%
intersections = Point2[]
for i=1:length(allFibres)-1
    for j=i+1:length(allFibres)
        testIntersection = intersects(allFibres[i],allFibres[j])
        if testIntersection[1]
            push!(intersections,testIntersection[2])
        end
    end
end


#%%

vertices = copy(intersections)
for f in allFibres
    append!(vertices,f.points)
end

#%%
for i in allFibres
    lines!(ax,i.points,color=:black)
end
scatter!(ax,intersections,color=:blue)
display(fig)