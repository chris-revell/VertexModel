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
includet("$(projectdir())/scripts/analysisFunctions/functions.jl")

dataDirectory = "data/old/2022-02-28-19-30-22"

# Import system data
conditionsDict    = load("$dataDirectory/dataFinal.jld2")
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$dataDirectory/matricesFinal.jld2")
@unpack A,B,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,externalF,boundaryVertices,cellPressures,edgeLengths = matricesDict["matrices"]

# Create vector of polygons for each cell
cellPolygons = makeCellPolygons(conditionsDict["params"],matricesDict["matrices"])

# Find cell midpoint links T
T = makeCellLinks(conditionsDict["params"],matricesDict["matrices"])

# Create vector of triangles from midpoint links
linkTriangles = makeLinkTriangles(conditionsDict["params"],matricesDict["matrices"])
linkTriangleAreas = abs.(area.(linkTriangles))

q = calculateSpokes(conditionsDict["params"],matricesDict["matrices"])

cellCurls = calculateCellCurls(conditionsDict["params"],matricesDict["matrices"])

cellDivs = calculateCellDivs(conditionsDict["params"],matricesDict["matrices"])

vertexDivs = calculateVertexDivs(conditionsDict["params"],matricesDict["matrices"],q,linkTriangleAreas)

vertexCurls = calculateVertexCurls(conditionsDict["params"],matricesDict["matrices"],q,linkTriangleAreas)

divLims = (-maximum(abs.([vertexDivs; cellDivs])),maximum(abs.([vertexDivs; cellDivs])))
curlLims = (-maximum(abs.(vertexCurls)),maximum(abs.(vertexCurls)))

# cellCurlLims = (-maximum(abs.(cellCurls)),maximum(abs.(cellCurls)))
# cellDivLims = (-maximum(abs.(cellDivs)),maximum(abs.(cellDivs)))
# vertexCurlLims = (-maximum(abs.(vertexCurls)),maximum(abs.(vertexCurls)))
# vertexDivLims = (-maximum(abs.(vertexDivs)),maximum(abs.(vertexDivs)))

# Set up figure canvas
set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica", fontsize=19)
fig = Figure(resolution=(1000,1500))
# grid = fig[1,1] = GridLayout()


# Cell div axis
ax1 = Axis(fig[1,1],aspect=DataAspect())
hidedecorations!(ax1)
hidespines!(ax1)
for i=1:nCells
    poly!(ax1,cellPolygons[i],color=[cellDivs[i]],colormap=:bwr,colorrange=divLims, strokecolor=(:black,1.0),strokewidth=1)
end
Label(fig[1,1,Bottom()],L"(a)",textsize = 32)
# Colorbar(fig[1,2][1,2],limits=cellDivLims,colormap=:bwr,flipaxis=true)

# Vertex div axis
ax2 = Axis(fig[1,2],aspect=DataAspect())
hidedecorations!(ax2)
hidespines!(ax2)
for k=1:nVerts
    poly!(ax2,linkTriangles[k],color=[vertexDivs[k]],colorrange=divLims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
end
for i=1:nCells
    poly!(ax2,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
end
Label(fig[1,2,Bottom()],L"(b)",textsize = 32)
Colorbar(fig[1,3],limits=divLims,colormap=:bwr,flipaxis=true)


# Cell curl axis
ax3 = Axis(fig[2,1],aspect=DataAspect())
hidedecorations!(ax3)
hidespines!(ax3)
for i=1:nCells
    poly!(ax3,cellPolygons[i],color=[cellCurls[i]],colormap=:bwr,colorrange=curlLims,strokecolor=(:black,1.0),strokewidth=1) #:bwr
end
Label(fig[2,1,Bottom()],L"(c)",textsize = 32)
# Colorbar(fig[1,1][1,2],limits=cellCurlLims,colormap=:bwr,flipaxis=true)

# Vertex curl axis
ax4 = Axis(fig[2,2],aspect=DataAspect())
hidedecorations!(ax4)
hidespines!(ax4)
for k=1:nVerts
    poly!(ax4,linkTriangles[k],color=[vertexCurls[k]],colorrange=curlLims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
end
for i=1:nCells
    poly!(ax4,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
end
Label(fig[2,2,Bottom()],L"(d)",textsize = 32)
Colorbar(fig[2,3],limits=curlLims,colormap=:bwr,flipaxis=true)



# Internal vertex curls
for k=1:nVerts
    boundaryVertices[k] != 0 ? vertexCurls[k]=0.0 : nothing
end
# Vertex couples
vertexCouples = Float64[]
for k=1:nVerts
    curl = 0
    if boundaryVertices[k] == 0
        for i=1:nCells
            for j=1:nEdges
                curl -= cellPressures[i]*B[i,j]*edgeLengths[j]^2*A[j,k]/(6.0*linkTriangleAreas[k])
            end
        end
    end
    push!(vertexCouples,curl)
end
otherLims = (-maximum(abs.([vertexCouples;vertexCurls])),maximum(abs.([vertexCouples;vertexCurls])))

# Internal vertex curl axis
ax5 = Axis(fig[3,1],aspect=DataAspect())
hidedecorations!(ax5)
hidespines!(ax5)
for k=1:nVerts
    poly!(ax5,linkTriangles[k],color=[vertexCurls[k]],colorrange=otherLims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
end
for i=1:nCells
    poly!(ax5,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
end
Label(fig[3,1,Bottom()],L"(e)",textsize = 32)

# Vertex couples axis
ax6 = Axis(fig[3,2],aspect=DataAspect())
hidedecorations!(ax6)
hidespines!(ax6)
for k=1:nVerts
    poly!(ax6,linkTriangles[k],color=[vertexCouples[k]],colorrange=otherLims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
end
for i=1:nCells
    poly!(ax6,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
end
Label(fig[3,2,Bottom()],L"(f)",textsize = 32)
cb = Colorbar(fig[3,3],limits=otherLims,colormap=:bwr,flipaxis=true) #:bwr
# cb.alignmode = Mixed(top=-10)

colsize!(fig.layout,1,Aspect(1,1.0))
colsize!(fig.layout,2,Aspect(1,1.0))


colgap!(fig.layout,1,Relative(0.0))
colgap!(fig.layout,2,Relative(0.0))
rowgap!(fig.layout,1,Relative(0.01))
rowgap!(fig.layout,2,Relative(0.01))


# boxes = [Box(g,color=(:white,0.0)) for g in [fig[1,1],fig[1,2],fig[1,3],fig[2,1],fig[2,2],fig[2,3],fig[3,1],fig[3,2],fig[3,3]]]

resize_to_layout!(fig)

display(fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/figure5.eps",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/figure5.png",fig)
