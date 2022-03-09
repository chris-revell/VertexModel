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
using LaTeXStrings

# Local modules
includet("$(projectdir())/src/VertexModelContainers.jl"); using .VertexModelContainers

function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end

dataDirectory = "data/sims/2022-02-28-19-30-22"

# Import system data
conditionsDict    = load("$dataDirectory/dataFinal.jld2")
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,pressureExternal,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
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

# Calculate div on each cell
cellDivs = Float64[]
for c=1:nCells
    cellVertices = findall(x->x!=0,C[c,:])
    vertexAngles = zeros(size(cellVertices))
    for (k,v) in enumerate(cellVertices)
        vertexAngles[k] = atan((R[v].-cellPositions[c])...)
    end
    m = minimum(vertexAngles)
    vertexAngles .-= m
    cellVertices .= cellVertices[sortperm(vertexAngles)]
    cellEdges = findall(x->x!=0,B[c,:])
    edgeAngles = zeros(size(cellEdges))
    for (k,e) in enumerate(cellEdges)
        edgeAngles[k] = atan((edgeMidpoints[e].-cellPositions[c])...)
    end
    edgeAngles .+= (2π-m)
    edgeAngles .= edgeAngles.%(2π)
    cellEdges .= cellEdges[sortperm(edgeAngles)]
    h = @SVector [0.0,0.0]
    divSum = 0
    for (i,e) in enumerate(cellEdges)
        h = h + ϵ*F[cellVertices[i],c]
        divSum -= B[c,e]*(h⋅(ϵ*edgeTangents[e]))/cellAreas[c]
    end
    divSum *= (-0.5)
    push!(cellDivs,divSum)
end

H = Diagonal(cellAreas)
boundaryEdgesFactor = abs.(boundaryEdges.-1)# =1 for internal vertices, =0 for boundary vertices
diagonalComponent = (boundaryEdgesFactor'.*((edgeLengths.^2)./(2.0.*trapeziumAreas)))[:,1] # Multiply by boundaryEdgesFactor vector to set boundary vertex contributions to zero
Tₑ = Diagonal(diagonalComponent)
Lf = (H\B)*Tₑ*Bᵀ
dropzeros!(Lf)

eigenvectors = (eigen(Matrix(Lf))).vectors
eigenvalues = (eigen(Matrix(Lf))).values

ḡ = ((onesVec'*H*cellDivs)/(onesVec'*H*ones(nCells))).*onesVec

ğ = cellDivs.-ḡ

ϕ̆ = zeros(nCells)
eigenmodeAmplitudes = Float64[]
for k=2:nCells
    numerator = eigenvectors[:,k]'*H*ğ
    denominator = eigenvalues[k]*(eigenvectors[:,k]'*H*eigenvectors[:,k])
    ϕ̆ .+= (numerator/denominator).*eigenvectors[:,k]
    push!(eigenmodeAmplitudes,(numerator/denominator))
end

ϕLims = (-maximum(abs.(ϕ̆)),maximum(abs.(ϕ̆)))

# update_theme!(font = "CMU Classical Serif", fontsize = 24)
fig = Figure(resolution=(1000,1000),fontsize = 24)
ax1 = Axis(fig[1,1][1,1],aspect=DataAspect(),title = L"$\phi^c$ Potential",fontsize=32)
hidedecorations!(ax1)
hidespines!(ax1)
# ax1.title = L"$\phi^c$ Potential"
for i=1:nCells
    poly!(ax1,cellPolygons[i],color=[ϕ̆[i]],colormap=:bwr,colorrange=ϕLims, strokecolor=(:black,1.0),strokewidth=5)
end
Colorbar(fig[1,1][1,2],limits=ϕLims,colormap=:bwr,flipaxis=false,align=:left)

ax2 = Axis(fig[2,:],title = LaTeXString("Eigenmode amplitudes"), xlabel=L"Eigenmode number, $i$", ylabel=LaTeXString("Amplitude"),fontsize=32)
lines!(ax2,collect(2:nCells),abs.(eigenmodeAmplitudes),linewidth=3)

display(fig)
# save("$dataDirectory/phi^cPotential.pdf",fig)
# save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/pdf/phi^cPotential.pdf",fig)
# save("$dataDirectory/phi^cPotential.svg",fig)
# save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/svg/phi^cPotential.svg",fig)
# save("$dataDirectory/phi^cPotential.png",fig)
# save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/png/phi^cPotential.png",fig)
#
