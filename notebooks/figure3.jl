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
using Printf

# Local modules
includet("$(projectdir())/scripts/analysisFunctions/functions.jl")

dataDirectory = "data/sims/2022-02-28-19-30-22"

# Import system data
conditionsDict    = load("$dataDirectory/dataFinal.jld2")
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,L₀,A₀,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$dataDirectory/matricesFinal.jld2")
@unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths = matricesDict["matrices"]

T = makeCellLinks(conditionsDict["params"],matricesDict["matrices"])

edgeTrapezia = makeEdgeTrapezia(conditionsDict["params"],matricesDict["matrices"])
trapeziumAreas = abs.(area.(edgeTrapezia))

linkTriangles = makeLinkTriangles(conditionsDict["params"],matricesDict["matrices"])
linkTriangleAreas = abs.(area.(linkTriangles))

cellPolygons = makeCellPolygons(conditionsDict["params"],matricesDict["matrices"])

Lf = makeLf(conditionsDict["params"],matricesDict["matrices"],trapeziumAreas)
decompositionLf = (eigen(Matrix(Lf))).vectors


Lᵥ = makeLv(conditionsDict["params"],matricesDict["matrices"],linkTriangleAreas,trapeziumAreas)
decompositionLv = (eigen(Matrix(Lᵥ))).vectors

# Set up figure canvas
set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
fig = Figure(resolution=(1500,1600))
superGrid = fig[1,1] = GridLayout()
grid1 = superGrid[1,1] = GridLayout()

for x=1:5
    for y=1:4
        eigenvectorIndex = ((y-1)*5 + x)+1
        lims = (-maximum(abs.(decompositionLv[:,eigenvectorIndex])),maximum(abs.(decompositionLv[:,eigenvectorIndex])))
        ax = Axis(grid1[y,x],aspect=DataAspect())
        hidedecorations!(ax)
        hidespines!(ax)
        for k=1:nVerts
            poly!(ax,linkTriangles[k],color=[decompositionLv[k,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
        end
        for i=1:nCells
            poly!(ax,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
        end
        Label(grid1[y,x,Bottom()],
                L"k=%$eigenvectorIndex",
                textsize = 16,
        )
    end
end
for x=1:5
    for y=1:4
        eigenvectorIndex = ((y-1)*5 + x)+(nVerts-20)
        lims = (-maximum(abs.(decompositionLv[:,eigenvectorIndex])),maximum(abs.(decompositionLv[:,eigenvectorIndex])))
        ax = Axis(grid1[y+4,x],aspect=DataAspect())
        hidedecorations!(ax)
        hidespines!(ax)
        for k=1:nVerts
            poly!(ax,linkTriangles[k],color=[decompositionLv[k,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
        end
        for i=1:nCells
            poly!(ax,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
        end
        Label(grid1[y+4,x,Bottom()],
                L"k=%$eigenvectorIndex",
                textsize = 16,
        )
    end
end
Label(grid1[9,:],
        L"(a)",
        textsize = 32,
)


vlineAx = Axis(superGrid[1,2])
hidedecorations!(vlineAx)
hidespines!(vlineAx)
vlines!(vlineAx,0.0,color=:black)
colsize!(superGrid, 2, Aspect(1, 0.01))


signInversions = [3, 5, 9, 10, 11, 12, 13, 17, 19, 20, 21]

grid2 = superGrid[1,3] = GridLayout()
for x=1:5
    for y=1:4
        eigenvectorIndex = ((y-1)*5 + x)+1
        lims = (-maximum(abs.(decompositionLf[:,eigenvectorIndex])),maximum(abs.(decompositionLf[:,eigenvectorIndex])))
        ax = Axis(grid2[y,x],aspect=DataAspect())
        hidedecorations!(ax)
        hidespines!(ax)
        if eigenvectorIndex in signInversions
            for i=1:nCells
                poly!(ax,cellPolygons[i],color=[-decompositionLf[i,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokecolor=(:black,1.0),strokewidth=1) #:bwr
            end
        else
            for i=1:nCells
                poly!(ax,cellPolygons[i],color=[decompositionLf[i,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokecolor=(:black,1.0),strokewidth=1) #:bwr
            end
        end
        Label(grid2[y,x,Bottom()],
                L"i=%$eigenvectorIndex",
                textsize = 16,
        )
    end
end
for x=1:5
    for y=1:4
        eigenvectorIndex = ((y-1)*5 + x)+(nCells-20)
        lims = (-maximum(abs.(decompositionLf[:,eigenvectorIndex])),maximum(abs.(decompositionLf[:,eigenvectorIndex])))
        ax = Axis(grid2[y+4,x],aspect=DataAspect())
        hidedecorations!(ax)
        hidespines!(ax)
        for i=1:nCells
            poly!(ax,cellPolygons[i],color=[decompositionLf[i,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokecolor=(:black,1.0),strokewidth=1) #:bwr
        end
        Label(grid2[y+4,x,Bottom()],
                L"i=%$eigenvectorIndex",
                textsize = 16,
        )
    end
end
Label(grid2[9,:],
        L"(b)",
        textsize = 32,
)

resize_to_layout!(fig)

# display(fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/figure3.eps",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/figure3.png",fig)
