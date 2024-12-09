using VertexModel
@from "$(srcdir("RotationMatrix.jl"))" using RotationMatrix
using LinearAlgebra
using GeometryBasics
using DrWatson 
using GLMakie
using Colors
using Random
using Distributions

using CairoMakie

R, matrices, params = loadData("24-11-27-14-47-29_L₀=0.75_nCells=61_realTimetMax=346000.0_γ=0.2",outputNumber=79)

@unpack nCells,
nEdges,
nVerts,
nonDimCycleTime,
distLogNormal,
γ, t1Threshold = params
@unpack A, 
B, 
cellTimeToDivide, 
cellPositions, 
edgeMidpoints, 
cellEdgeCount, 
cellPerimeters,
cellVertexOrders, 
cellEdgeOrders,
cellPerpAxes, 
boundaryEdges, 
edgeTangents,
μ, 
Γ = matrices


errorCells = Int64[]
shortAxisLines = []


fig = CairoMakie.Figure()
ax = Axis(fig[1,1], aspect=DataAspect()) 

for i=1:nCells
    spokes = [R[kk].-matrices.cellPositions[i] for kk in matrices.cellVertexOrders[i][0:end]]
    
    crossVec = matrices.cellPerpAxes[i]×[1,0,0]
    ϵCoordinates = ϵ(v=crossVec, θ=asin(norm(crossVec)/(norm(matrices.cellPerpAxes[i]))))
    rotatedSpokes = [(ϵCoordinates*s)[2:end] for s in spokes]
    cellShapeTensor = sum(rotatedSpokes[2:end].*transpose.(rotatedSpokes[2:end]))./matrices.cellEdgeCount[i]
    
    # Long and short axis from eigenvectors of shapetensor
    # Put some sort of tolerance that if eigenvalues are approx equal we randomly choose a division orientation, eg circ >0.95
    eigenVals, eigenVecs = eigen(cellShapeTensor) # eigenvalues/vectors listed smallest to largest eigval.
    circ = abs(eigenVals[1]/eigenVals[2]) # circularity
    
    # #for very circular cells randomly choose division axis
    # if circ > 0.95 
    #     theta = rand()*π
    #     shortvec = [cos(theta), sin(theta)]
    # else
        
    if eigenVecs[:,1][2] < 0.0 #make it so vector is pointing in positive y direction (to fit with existing code in assigning new edges)
        shortvec = -1.0*cellPerimeters[i].*eigenVecs[:,1]
    else
        shortvec = cellPerimeters[i].*eigenVecs[:,1]
    end

    # end
    shortAxisLine = Line(Point{2,Float64}(shortvec), Point{2,Float64}(-shortvec))
    push!(shortAxisLines, shortvec)
    # Test cell edges for an intersection
    poly = LineString(Point{2, Float64}.(rotatedSpokes)) # Start and end with the same vertex by indexing circular array from 0 to end
    intersections = [intersects(line, shortAxisLine) for line in poly] #find which edges intersect and where
    intersectedIndices = findall(x->x!=0, first.(intersections))
    if length(intersectedIndices) < 2 
        push!(errorCells,i)
        @show circ
        poly!(ax, first.(rotatedSpokes[1:end-1]), last.(rotatedSpokes[1:end-1]))
        lines!(ax, first.([shortvec, -shortvec]), last.([shortvec, -shortvec]))
        display(fig)
        sleep(2)
        empty!(ax)
    end

end

display(fig)










    # set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    fig = GLMakie.Figure(size=(1000,1000))
    ax = LScene(fig[1,1])

    for i in errorCells #i=1:nCells
        verts = Float64[]
        for k=1:length(cellVertexOrders[i])
            append!(verts, R[cellVertexOrders[i][k]])
            append!(verts, R[cellVertexOrders[i][k+1]])
            append!(verts, cellPositions[i])
        end
        connectedVerts = connect(verts, Point{3})
        connectedFaces = connect(1:length(connectedVerts), TriangleFace)
        mesh!(ax, connectedVerts, connectedFaces, color=RGB(rand(Xoshiro(i),3)...), shading=NoShading)
    end 

    for i in errorCells
        lines!(ax, Point3.([cellPositions[i].+shortAxisLines[i], cellPositions[i].-shortAxisLines[i]]))
    end

    if labels
        text!(ax, Point{3,Float64}.([p.+[0.0,0.0,0.1] for p in matrices.cellPositions]), text = string.(collect(1:nCells)), color=:red)
        scatter!(ax, Point{3,Float64}.([p.+[0.0,0.0,0.1] for p in matrices.cellPositions]), color=:red)
        text!(ax, Point{3,Float64}.([p.+[0.0,0.0,0.1] for p in matrices.edgeMidpoints]), text = string.(collect(1:nEdges)), color=:green)
        scatter!(ax, Point{3,Float64}.([p.+[0.0,0.0,0.1] for p in matrices.edgeMidpoints]), color=:green)
        text!(ax, Point{3,Float64}.([p.+[0.0,0.0,0.1] for p in R]), text = string.(collect(1:nVerts)), color=:blue)
        scatter!(ax, Point{3,Float64}.([p.+[0.0,0.0,0.1] for p in R]), color=:blue)
    end
    reset_limits!(ax)
    
    display(fig)