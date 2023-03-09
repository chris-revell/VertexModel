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

function allCurlsDivs(dataDirectory, show)
    isdir("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/png") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/png")
    isdir("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf")
    isdir("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg")
    isdir("$dataDirectory/png") ? nothing : mkpath("$dataDirectory/png")
    isdir("$dataDirectory/pdf") ? nothing : mkpath("$dataDirectory/pdf")
    isdir("$dataDirectory/svg") ? nothing : mkpath("$dataDirectory/svg")

    # Import system data
    conditionsDict    = load("$dataDirectory/dataFinal.jld2")
    @unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
    matricesDict = load("$dataDirectory/matricesFinal.jld2")
    @unpack A,B,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,externalF,boundaryVertices = matricesDict["matrices"]

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

    # Set up figure canvas
    fig = Figure(resolution=(900,1000))
    grid = fig[1,1] = GridLayout()

    # Cell div axis
    ax1 = Axis(grid[1,1][1,1],aspect=DataAspect())
    hidedecorations!(ax1)
    hidespines!(ax1)
    for i=1:nCells
        poly!(ax1,cellPolygons[i],color=[cellDivs[i]],colormap=:bwr,colorrange=divLims, strokecolor=(:black,1.0),strokewidth=5)
    end
    Label(grid[1, 1, Bottom()],LaTeXString("(a)"),fontsize = 28)
    # Label(grid[1, 1, Left()],LaTeXString("div"),fontsize = 28)
    # Label(grid[1, 1, Top()],LaTeXString("Cells"),fontsize = 28)

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
    Label(grid[1, 2, Bottom()],LaTeXString("(b)"),fontsize = 28)
    # Label(grid[1, 2, Top()],LaTeXString("Vertices"),fontsize = 28)
    Colorbar(grid[1,3],limits=divLims,colormap=:bwr,flipaxis=false)

    # Cell curl axis
    ax3 = Axis(grid[2,1],aspect=DataAspect())
    hidedecorations!(ax3)
    hidespines!(ax3)
    for i=1:nCells
        poly!(ax3,cellPolygons[i],color=[cellCurls[i]],colormap=:bwr,colorrange=curlLims,strokecolor=(:black,1.0),strokewidth=5) #:bwr
    end
    Label(grid[2, 1, Bottom()],LaTeXString("(c)"),fontsize = 28)
    # Label(grid[2, 1, Left()],LaTeXString("curl"),fontsize = 28)

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
    Label(grid[2, 2, Bottom()],LaTeXString("(d)"),fontsize = 28)
    Colorbar(grid[2,3],limits=curlLims,colormap=:bwr,flipaxis=false)

    show==1 ? display(fig) : nothing
    save("$dataDirectory/pdf/allCurlsDivs.pdf",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/allCurlsDivs.pdf",fig)
    save("$dataDirectory/svg/allCurlsDivs.svg",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/allCurlsDivs.svg",fig)
    save("$dataDirectory/png/allCurlsDivs.png",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/allCurlsDivs.png",fig)
end
