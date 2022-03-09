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
@unpack A,B,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,externalF = matricesDict["matrices"]

# Create vector of polygons for each cell
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

# Find cell midpoint links T
onesVec = ones(1,nCells)
boundaryEdges = abs.(onesVec*B)
cᵖ = boundaryEdges'.*edgeMidpoints
T = SVector{2,Float64}[]
for j=1:nEdges
    Tⱼ = @SVector zeros(2)
    for i=1:nCells
        Tⱼ = Tⱼ + B[i,j]*(cellPositions[i].-cᵖ[j])
    end
    push!(T,Tⱼ)
end

# Create vector of triangles from midpoint links
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

# Calculate curl on each cell
cellCurls = Float64[]
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
    curlSum = 0
    for (i,e) in enumerate(cellEdges)
        h = h + ϵ*F[cellVertices[i],c]
        curlSum += B[c,e]*(h⋅edgeTangents[e])/cellAreas[c]
    end
    push!(cellCurls,curlSum)
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
        # pressureAngle = atan((-1.0.*externalF[k])...)
        # vertexEdges = findall(x->x!=0,A[:,k])
        # edgeAngles = zeros(length(vertexEdges))
        # for (i,e) in enumerate(vertexEdges)
        #     edgeAngles[i] = atan((-A[e,k].*edgeTangents[e])...)
        # end
        #
        # edgeAngles .+= 2π-pressureAngle
        # edgeAngles .= edgeAngles.%2π
        # vertexEdges .= vertexEdges[sortperm(edgeAngles,rev=true)]
        #
        # vertexCells = findall(x->x!=0,C[:,k])
        # cellAngles = zeros(length(vertexCells))
        # for i=1:length(cellAngles)
        #     cellAngles[i] = atan((cellPositions[vertexCells[i]].-R[k])...)
        # end
        # cellAngles .+= 2π-pressureAngle
        # cellAngles .= cellAngles.%2π
        # vertexCells .= vertexCells[sortperm(cellAngles,rev=true)]
        #
        # h = @SVector [0.0,0.0]
        divSum = 0
        # h = h + ϵ*externalF[k]
        # divSum -= A[vertexEdges[1],k]*((ϵₖ*T[vertexEdges[1]])⋅h)/E[k]
        # for (i,j) in enumerate(vertexEdges[2:end])
        #     h = h + ϵ*F[k,vertexCells[i]]
        #     divSum -= A[j,k]*((ϵₖ*T[j])⋅h)/E[k]
        # end
    end
    divSum *= (-0.5)
    push!(vertexDivs,divSum)
end



# Calculate curl at each vertex
vertexCurls = Float64[]
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
        curlSum = 0
        for (i,j) in enumerate(vertexEdges)
            h = h + ϵ*F[k,vertexCells[i]]
            curlSum += A[j,k]*(T[j]⋅h)/E[k]
        end
    else
        # # Set angles relative to pressure force angle, equivalent to the angle of a cell that doesn't actually exist
        # pressureAngle = atan((-1.0.*externalF[k])...)
        # vertexEdges = findall(x->x!=0,A[:,k])
        # edgeAngles = zeros(length(vertexEdges))
        # for (i,e) in enumerate(vertexEdges)
        #     edgeAngles[i] = atan((-A[e,k].*edgeTangents[e])...)
        # end
        #
        # edgeAngles .+= 2π-pressureAngle
        # edgeAngles .= edgeAngles.%2π
        # vertexEdges .= vertexEdges[sortperm(edgeAngles,rev=true)]
        #
        # vertexCells = findall(x->x!=0,C[:,k])
        # cellAngles = zeros(length(vertexCells))
        # for i=1:length(cellAngles)
        #     cellAngles[i] = atan((cellPositions[vertexCells[i]].-R[k])...)
        # end
        # cellAngles .+= 2π-pressureAngle
        # cellAngles .= cellAngles.%2π
        # vertexCells .= vertexCells[sortperm(cellAngles,rev=true)]
        #
        # h = @SVector [0.0,0.0]
        # curlSum = 0
        # h = h + ϵ*externalF[k]
        # curlSum += A[vertexEdges[1],k]*(T[vertexEdges[1]]⋅h)/E[k]
        # for (i,j) in enumerate(vertexEdges[2:end])
        #     h = h + ϵ*F[k,vertexCells[i]]
        #     curlSum += A[j,k]*(T[j]⋅h)/E[k]
        # end
        curlSum = 0.0
    end
    push!(vertexCurls,curlSum)
end

divLims = (-maximum(abs.([vertexDivs; cellDivs])),maximum(abs.([vertexDivs; cellDivs])))
curlLims = (-maximum(abs.(vertexCurls)),maximum(abs.(vertexCurls)))

# Set up figure canvas
fig = Figure(resolution=(900,1000))
grid = fig[1,1] = GridLayout()
# Cell div axis
ax1 = Axis(grid[1,1][1,1],aspect=DataAspect())
hidedecorations!(ax1)
hidespines!(ax1)
# ax1.title = "-0.5 x Cell divs"
for i=1:nCells
    poly!(ax1,cellPolygons[i],color=[cellDivs[i]],colormap=:bwr,colorrange=divLims, strokecolor=(:black,1.0),strokewidth=5)
end
Label(grid[1, 1, Bottom()],
        LaTeXString("(a)"),
        textsize = 34,
        #halign = :right
)
Label(grid[1, 1, Left()],
        LaTeXString("div"),
        textsize = 34,
        #halign = :right
)
Label(grid[1, 1, Top()],
        LaTeXString("Cells"),
        textsize = 34,
        #halign = :right
)

# Vertex div axis
ax2 = Axis(grid[1,2],aspect=DataAspect())
hidedecorations!(ax2)
hidespines!(ax2)
for k=1:nVerts
    poly!(ax2,linkTriangles[k],color=[vertexDivs[k]],colorrange=divLims,colormap=:bwr,strokewidth=2,strokecolor=(:black,0.0)) #:bwr
end
for i=1:nCells
    poly!(ax2,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=5) #:bwr
end
Label(grid[1, 2, Bottom()],
        LaTeXString("(b)"),
        textsize = 34,
        #halign = :right
)
Label(grid[1, 2, Top()],
        LaTeXString("Vertices"),
        textsize = 34,
        #halign = :right
)
Colorbar(grid[1,3],limits=divLims,colormap=:bwr,flipaxis=false)

# Cell curl axis
ax3 = Axis(grid[2,1],aspect=DataAspect())
hidedecorations!(ax3)
hidespines!(ax3)
for i=1:nCells
    poly!(ax3,cellPolygons[i],color=[cellCurls[i]],colormap=:bwr,colorrange=curlLims,strokecolor=(:black,1.0),strokewidth=5) #:bwr
end
Label(grid[2, 1, Bottom()],
        LaTeXString("(c)"),
        textsize = 34,
)
Label(grid[2, 1, Left()],
        LaTeXString("curl"),
        textsize = 34,
)

# Vertex curl axis
ax4 = Axis(grid[2,2][1,1],aspect=DataAspect())
hidedecorations!(ax4)
hidespines!(ax4)
for k=1:nVerts
    poly!(ax4,linkTriangles[k],color=[vertexCurls[k]],colorrange=curlLims,colormap=:bwr,strokewidth=2,strokecolor=(:black,0.0)) #:bwr
end
for i=1:nCells
    poly!(ax4,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=5) #:bwr
end
Label(grid[2, 2, Bottom()],
        LaTeXString("(d)"),
        textsize = 34,
)
Colorbar(grid[2,3],limits=curlLims,colormap=:bwr,flipaxis=false)


# colgap!(grid,-10)

display(fig)
save("$dataDirectory/allCurlsDivs.pdf",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/pdf/allCurlsDivs.pdf",fig)
save("$dataDirectory/allCurlsDivs.svg",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/svg/allCurlsDivs.svg",fig)
save("$dataDirectory/allCurlsDivs.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/png/allCurlsDivs.png",fig)
