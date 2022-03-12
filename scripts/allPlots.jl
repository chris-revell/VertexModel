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
using Printf

# Local modules
includet("$(projectdir())/src/VertexModelContainers.jl"); using .VertexModelContainers

function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end

dataDirectory = "data/sims/2022-03-10-21-07-10"

# Import system data
conditionsDict    = load("$dataDirectory/dataFinal.jld2")
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,pressureExternal,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$dataDirectory/matricesFinal.jld2")
@unpack A,Aᵀ,B,Bᵀ,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,externalF,boundaryVertices,edgeLengths = matricesDict["matrices"]

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
            divSum -= A[j,k]*((ϵₖ*T[j])⋅h)/linkTriangleAreas[k]
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
            curlSum += A[j,k]*(T[j]⋅h)/linkTriangleAreas[k]
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

# Laplacians
H = Diagonal(cellAreas)
boundaryEdgesFactor = abs.(boundaryEdges.-1)# =1 for internal vertices, =0 for boundary vertices
diagonalComponent = (boundaryEdgesFactor'.*((edgeLengths.^2)./(2.0.*trapeziumAreas)))[:,1] # Multiply by boundaryEdgesFactor vector to set boundary vertex contributions to zero
Tₑ = Diagonal(diagonalComponent)
Lf = (H\B)*Tₑ*Bᵀ
dropzeros!(Lf)
E = Diagonal(linkTriangleAreas)
Tₑ = Diagonal((edgeLengths.^2)./(2.0.*trapeziumAreas))
Lᵥ = (E\Aᵀ)*(Tₑ\A)
dropzeros!(Lᵥ)



# allCurlsDivs figure
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
save("$dataDirectory/allCurlsDivs.pdf",fig)
#save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/allCurlsDivs.pdf",fig)
save("$dataDirectory/allCurlsDivs.svg",fig)
#save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/allCurlsDivs.svg",fig)
save("$dataDirectory/allCurlsDivs.png",fig)
#save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/allCurlsDivs.png",fig)

# All eigenmodes Lf
decomposition = (eigen(Matrix(Lf))).vectors
isdir("$dataDirectory/eigenmodesLf") ? nothing : mkpath("$dataDirectory/eigenmodesLf")
isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/eigenmodesLf") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/eigenmodesLf")
# Set up figure canvas
fig = Figure(resolution=(1000,1000))
grid = fig[1,1] = GridLayout()
ax = Axis(grid[1,1],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
for eigenvectorIndex=2:nCells
    empty!(ax)
    lims = (minimum(decomposition[:,eigenvectorIndex]),maximum(decomposition[:,eigenvectorIndex]))
    for i=1:nCells
        poly!(ax,cellPolygons[i],color=[decomposition[i,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokecolor=(:black,1.0),strokewidth=1) #:bwr
    end
    save("$dataDirectory/eigenmodesLf/eigenmode$(@sprintf("%03d", eigenvectorIndex)).png",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/eigenmodesLf/eigenmode$(@sprintf("%03d", eigenvectorIndex)).png",fig)
end

# All eigenmodes Lv
decomposition = (eigen(Matrix(Lᵥ))).vectors
isdir("$dataDirectory/eigenmodesLv") ? nothing : mkpath("$dataDirectory/eigenmodesLv")
isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/eigenmodesLv") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/eigenmodesLv")
fig = Figure(resolution=(1000,1000))
grid = fig[1,1] = GridLayout()
ax = Axis(grid[1,1],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
for eigenvectorIndex=2:nVerts
    empty!(ax)
    lims = (minimum(decomposition[:,eigenvectorIndex]),maximum(decomposition[:,eigenvectorIndex]))
    for k=1:nVerts
        poly!(ax,linkTriangles[k],color=[decomposition[k,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
    end
    for i=1:nCells
        poly!(ax,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
    end
    save("$dataDirectory/eigenmodesLv/eigenmode$(@sprintf("%03d", eigenvectorIndex)).png",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/eigenmodesLv/eigenmode$(@sprintf("%03d", eigenvectorIndex)).png",fig)
end

# Force neighbourhood
centralCell=14
# Find all cells neighbouring original cell
cellNeighbourMatrix = B*Bᵀ
dropzeros!(cellNeighbourMatrix)
neighbouringCells = findall(!iszero,cellNeighbourMatrix[centralCell,:])
# Find and sort all vertices around cells neighbouring centralCell
cellVerticesDict = Dict()
for c in neighbouringCells
    # Find vertices around cell
    cellVertices = findall(x->x!=0,C[c,:])
    # Find angles of vertices around cell
    vertexAngles = zeros(size(cellVertices))
    for (k,v) in enumerate(cellVertices)
        vertexAngles[k] = atan((R[v].-cellPositions[c])...)
    end
    # Sort vertices around cell by polar angle
    cellVertices .= cellVertices[sortperm(vertexAngles)]
    # Store sorted cell vertices for this cell
    cellVerticesDict[c] = cellVertices
end
edgeMidpointPolygons = Vector{Point2f}[]
for i=1:nCells
    cellEdges = findall(!iszero,B[i,:])
    edgeAngles = zeros(size(cellEdges))
    for (k,v) in enumerate(cellEdges)
        edgeAngles[k] = atan((edgeMidpoints[v].-cellPositions[i])...)
    end
    cellEdges .= cellEdges[sortperm(edgeAngles)]
    push!(edgeMidpointPolygons,Point2f.(edgeMidpoints[cellEdges]))
end
fig = Figure(resolution=(1000,1000))
ax1 = Axis(fig[1,1],aspect=DataAspect())
hidedecorations!(ax1)
hidespines!(ax1)
for c in neighbouringCells
    poly!(ax1,cellPolygons[c],color=(getRandomColor(c),0.50),strokecolor=:black,strokewidth=8)
    poly!(ax1,edgeMidpointPolygons[c],color=(:white,0.0),strokecolor=(getRandomColor(c),1.0),strokewidth=8)
    scatter!(ax1,Point2f.(R[cellVerticesDict[c]]),color=:blue,markersize=16)
    scatter!(ax1,Point2f.(edgeMidpoints[findall(!iszero,B[c,:])]),color=:green,markersize=16)
end
scatter!(ax1,Point2f.(cellPositions[push!(neighbouringCells,centralCell)]),color=:red,markersize=16)
save("$dataDirectory/cell$(centralCell)ForceNeighbourhood.pdf",fig)
#save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/cell$(centralCell)ForceNeighbourhood.pdf",fig)
save("$dataDirectory/cell$(centralCell)ForceNeighbourhood.svg",fig)
#save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/cell$(centralCell)ForceNeighbourhood.svg",fig)
save("$dataDirectory/cell$(centralCell)ForceNeighbourhood.png",fig)
#save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/cell$(centralCell)ForceNeighbourhood.png",fig)

# Force network
fig = Figure(resolution=(1000,1000))
ax2 = Axis(fig[1,1],aspect=DataAspect())
hidedecorations!(ax2)
hidespines!(ax2)
centralCellVertices = findall(x->x!=0,C[centralCell,:])
centralVertexAngles = zeros(size(centralCellVertices))
for (k,v) in enumerate(centralCellVertices)
   centralVertexAngles[k] = atan((R[v].-cellPositions[centralCell])...)
end
m = minimum(centralVertexAngles)
centralVertexAngles .-= m
centralCellVertices .= centralCellVertices[sortperm(centralVertexAngles)]
# Sort cells neighbouring centralCell by angle
setdiff!(neighbouringCells,[centralCell]) # Remove centralCell from neighbours list
neighbourAngles = zeros(length(neighbouringCells))
for (i,c) in enumerate(neighbouringCells)
    neighbourAngles[i] = atan((cellPositions[c].-cellPositions[centralCell])...)
end
neighbourAngles .+= (2π-m)
neighbourAngles = neighbourAngles.%(2π)
neighbouringCells .= neighbouringCells[sortperm(neighbourAngles)]
# Draw force network
startPosition = @SVector [0.0,0.0]
for (i,v) in enumerate(cellVerticesDict[centralCell])
    arrows!(ax2,Point2f.([startPosition]),Vec2f.([ϵ*F[v,centralCell]]),linewidth=6,arrowsize=30,color=(getRandomColor(centralCell),0.9))
    #annotations!(ax2,string.([v]),Point2f.([startPosition.+ϵ*F[v,centralCell]./2.0]),color=(getRandomColor(centralCell),0.75))
    startPosition = startPosition + ϵ*F[v,centralCell]
    H = Array{SVector{2,Float64}}(undef,length(cellVerticesDict[neighbouringCells[i]])+1)
    cellForces = SVector{2, Float64}[]
    # Circular permutation of vertices to ensure vertex v is the first index
    # in the ordered cellVertices list around cell neighbouringCells[i]
    index = findall(x->x==v, cellVerticesDict[neighbouringCells[i]])
    cellVertices = circshift(cellVerticesDict[neighbouringCells[i]],1-index[1])
    H[1] = startPosition
    for (j,cv) in enumerate(cellVertices)
        push!(cellForces,+ϵ*F[cv,neighbouringCells[i]])
        H[j+1] = H[j]+cellForces[end]
    end
    arrows!(ax2,Point2f.(H),Vec2f.(cellForces),color=(getRandomColor(neighbouringCells[i]),0.9),linewidth=6,arrowsize=30)
end
save("$dataDirectory/cell$(centralCell)ForceNetwork.pdf",fig)
#save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/cell$(centralCell)ForceNetwork.pdf",fig)
save("$dataDirectory/cell$(centralCell)ForceNetwork.svg",fig)
#save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/cell$(centralCell)ForceNetwork.svg",fig)
save("$dataDirectory/cell$(centralCell)ForceNetwork.png",fig)
#save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/cell$(centralCell)ForceNetwork.png",fig)

# Full system
fig = Figure(resolution=(1000,1000))
ax = Axis(fig[1,1],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
for i=1:nCells
    if cellNeighbourMatrix[centralCell,i] == 0
        poly!(ax,cellPolygons[i],color=(getRandomColor(i),0.25),strokecolor=(:black,0.5),strokewidth=1)
    else
        poly!(ax,cellPolygons[i],color=(getRandomColor(i),1.0),strokecolor=(:black,1.0),strokewidth=2)
    end
end
for j=1:nEdges
    edgeCells = findall(!iszero,B[:,j])
    if boundaryEdges[j] == 0
        lines!(ax,Point2f.(cellPositions[edgeCells]),linewidth=1,color=(:white,1.0))
    else
        lines!(ax,Point2f.([cellPositions[edgeCells[1]],cᵖ[j]]),linewidth=1,color=(:white,1.0))
    end
end
scatter!(ax,Point2f.(R),alpha=0.5,color=:blue)
scatter!(ax,Point2f.(cellPositions),color=:red)
save("$dataDirectory/fullSystem.pdf",fig)
#save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/fullSystem.pdf",fig)
save("$dataDirectory/fullSystem.svg",fig)
#save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/fullSystem.svg",fig)
save("$dataDirectory/fullSystem.png",fig)
#save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/fullSystem.png",fig)
