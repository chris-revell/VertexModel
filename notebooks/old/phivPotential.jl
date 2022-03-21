# Import Julia packages
using DrWatson
@quickactivate
using Revise
using LinearAlgebra
using DelimitedFiles
using SparseArrays
using StaticArrays
using CairoMakie
using UnPack
using GeometryBasics
using Random
using Colors
using JLD2

# Local modules
includet("$(projectdir())/src/VertexModelContainers.jl"); using .VertexModelContainers

function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end

dataDirectory = "data/sims/2022-02-28-19-30-22"

# Import system data
conditionsDict    = load("$dataDirectory/dataFinal.jld2")
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$dataDirectory/matricesFinal.jld2")
@unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths,externalF = matricesDict["matrices"]

onesVec = ones(nCells)
boundaryEdges = abs.(onesVec'*B)
cᵖ = boundaryEdges'.*edgeMidpoints
T = SVector{2,Float64}[]
for j=1:nEdges
    Tⱼ = @SVector zeros(2)
    for i=1:nCells
        Tⱼ = Tⱼ + B[i,j]*(cellPositions[i].-cᵖ[j])
    end
    push!(T,Tⱼ)
end

edgeTrapezia = Vector{Point2f}[]
for j=1:nEdges
    edgeCells = findall(x->x!=0,B[:,j])
    edgeVertices = findall(x->x!=0,A[j,:])
    trapeziumVertices = [R[edgeVertices]; cellPositions[edgeCells]]
    com = sum(trapeziumVertices)./length(trapeziumVertices)
    angles = Float64[]
    for p=1:length(trapeziumVertices)
        angle = atan((trapeziumVertices[p].-com)...)
        push!(angles,angle)
    end
    trapeziumVertices .= trapeziumVertices[sortperm(angles)]
    push!(edgeTrapezia,Point2f.(trapeziumVertices))
end
trapeziumAreas = abs.(area.(edgeTrapezia))

linkTriangles = Vector{Point2f}[]
for k=1:nVerts
    if boundaryVertices[k] == 0
        vertexCells = findall(x->x!=0,C[:,k])
        push!(linkTriangles, Point2f.(cellPositions[vertexCells]))
    else
        vertexCells = findall(x->x!=0,C[:,k])
        vertexEdges = findall(x->x!=0,A[:,k])
        boundaryVertexEdges = intersect(vertexEdges,findall(x->x!=0,boundaryEdges[1,:]))
        kiteVertices = [edgeMidpoints[boundaryVertexEdges]; cellPositions[vertexCells]]
        push!(kiteVertices,R[k])
        com = sum(kiteVertices)./length(kiteVertices)
        angles = Float64[]
        for p=1:length(kiteVertices)
            angle = atan((kiteVertices[p].-com)...)
            push!(angles,angle)
        end
        kiteVertices .= kiteVertices[sortperm(angles)]
        push!(linkTriangles,Point2f.(kiteVertices))
    end
end
linkTriangleAreas = abs.(area.(linkTriangles))
E = linkTriangleAreas

cellPolygons = Vector{Point2f}[]
for i=1:nCells
    cellVertices = findall(x->x!=0,C[i,:])
    vertexAngles = zeros(size(cellVertices))
    for (k,v) in enumerate(cellVertices)
        vertexAngles[k] = atan((R[v].-cellPositions[i])...)
    end
    cellVertices .= cellVertices[sortperm(vertexAngles)]
    push!(cellPolygons,Point2f.(R[cellVertices]))
end

# Rotation matrix around vertices is the opposite of that around cells
ϵₖ = -1*ϵ

# Calculate div at each vertex
vertexDivs = Float64[]
# Working around a given vertex, an h force space point from a cell is mapped to the next edge anticlockwise from the cell
for k=1:nVerts
    if boundaryVertices[k] == 0
        vertexEdges = findall(x->x!=0,A[:,k])
        edgeAngles = zeros(length(vertexEdges))
        for (i,e) in enumerate(vertexEdges)
            edgeAngles[i] = atan((-A[e,k].*edgeTangents[e])...)
        end
        m = minimum(edgeAngles)
        edgeAngles .-= m
        vertexEdges .= vertexEdges[sortperm(edgeAngles,rev=true)]

        vertexCells = findall(x->x!=0,C[:,k])
        cellAngles = zeros(length(vertexCells))
        for i=1:length(cellAngles)
            cellAngles[i] = atan((cellPositions[vertexCells[i]].-R[k])...)
        end
        cellAngles .+= 2π-m
        cellAngles .= cellAngles.%2π
        vertexCells .= vertexCells[sortperm(cellAngles,rev=true)]
        h = @SVector [0.0,0.0]
        divSum = 0
        for (i,j) in enumerate(vertexEdges)
            h = h + ϵ*F[k,vertexCells[i]]
            divSum -= A[j,k]*((ϵₖ*T[j])⋅h)/E[k]
        end
    else
        # Set angles relative to pressure force angle, equivalent to the angle of a cell that doesn't actually exist
        pressureAngle = atan((-1.0.*externalF[k])...)
        vertexEdges = findall(x->x!=0,A[:,k])
        edgeAngles = zeros(length(vertexEdges))
        for (i,e) in enumerate(vertexEdges)
            edgeAngles[i] = atan((-A[e,k].*edgeTangents[e])...)
        end

        edgeAngles .+= 2π-pressureAngle
        edgeAngles .= edgeAngles.%2π
        vertexEdges .= vertexEdges[sortperm(edgeAngles,rev=true)]

        vertexCells = findall(x->x!=0,C[:,k])
        cellAngles = zeros(length(vertexCells))
        for i=1:length(cellAngles)
            cellAngles[i] = atan((cellPositions[vertexCells[i]].-R[k])...)
        end
        cellAngles .+= 2π-pressureAngle
        cellAngles .= cellAngles.%2π
        vertexCells .= vertexCells[sortperm(cellAngles,rev=true)]

        h = @SVector [0.0,0.0]
        divSum = 0
        h = h + ϵ*externalF[k]
        divSum -= A[vertexEdges[1],k]*((ϵₖ*T[vertexEdges[1]])⋅h)/E[k]
        for (i,j) in enumerate(vertexEdges[2:end])
            h = h + ϵ*F[k,vertexCells[i]]
            divSum -= A[j,k]*((ϵₖ*T[j])⋅h)/E[k]
        end
    end
    divSum *= (-0.5)
    push!(vertexDivs,divSum)
end

E = Diagonal(linkTriangleAreas)
Tₑ = Diagonal((edgeLengths.^2)./(2.0.*trapeziumAreas))
Lᵥ = (E\Aᵀ)*(Tₑ\A)
dropzeros!(Lᵥ)

eigenvectors = (eigen(Matrix(Lᵥ))).vectors
eigenvalues = (eigen(Matrix(Lᵥ))).values

ḡ = ((onesVec'*E*vertexDivs)/(onesVec'*E*ones(nCells))).*onesVec

ğ = cellDivs.-ḡ

ϕ̆ = zeros(nVerts)
eigenmodeAmplitudes = Float64[]
for k=2:nVerts
    numerator = eigenvectors[:,k]'*H*ğ
    denominator = eigenvalues[k]*(eigenvectors[:,k]'*E*eigenvectors[:,k])
    ϕ̆ .+= (numerator/denominator).*eigenvectors[:,k]
    push!(eigenmodeAmplitudes,(numerator/denominator))
end

ϕLims = (-maximum(abs.(ϕ̆)),maximum(abs.(ϕ̆)))

fig = Figure(resolution=(1000,1000))
grid = fig[1,1] = GridLayout()
ax1 = Axis(grid[1,1],aspect=DataAspect())
hidedecorations!(ax1)
hidespines!(ax1)
ax1.title = "ϕᵛ Potential"
for k=1:nVerts
    poly!(ax1,linkTriangles[k],color=[ϕ̆[k]],colorrange=ϕLims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
end
for i=1:nCells
    poly!(ax1,cellPolygons[i],color=(:white,0.0), strokecolor=(:black,1.0),strokewidth=5)
end
Colorbar(grid[1,2],limits=ϕLims,colormap=:bwr,flipaxis=false)

ax2 = Axis(grid[2,:])
lines!(ax2,collect(2:nVerts),abs.(eigenmodeAmplitudes),linewidth=3)
ax2.title = "Eigenmode amplitudes"
ax2.xlabel = "i"
ax2.ylabel = "Amplitude"

display(fig)
save("$dataDirectory/phi^vPotential.pdf",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/phi^vPotential.pdf",fig)
save("$dataDirectory/phi^vPotential.svg",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/phi^vPotential.svg",fig)
save("$dataDirectory/phi^vPotential.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/phi^vPotential.png",fig)
