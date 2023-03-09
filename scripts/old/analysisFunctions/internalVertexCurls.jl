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

function internalVertexCurls(dataDirectory, show)
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

    vertexCurls = calculateVertexCurls(conditionsDict["params"],matricesDict["matrices"],q,linkTriangleAreas)

    for k=1:nVerts
        boundaryVertices[k] != 0 ? vertexCurls[k]=0.0 : nothing
    end

    curlLims = (-maximum(abs.(vertexCurls)),maximum(abs.(vertexCurls)))

    # Set up figure canvas
    fig = Figure(resolution=(900,1000))
    grid = fig[1,1] = GridLayout()

    # Vertex curl axis
    ax4 = Axis(grid[1,1],aspect=DataAspect())
    hidedecorations!(ax4)
    hidespines!(ax4)
    for k=1:nVerts
        poly!(ax4,linkTriangles[k],color=[vertexCurls[k]],colorrange=curlLims,colormap=:bwr,strokewidth=2,strokecolor=(:black,0.0)) #:bwr
    end
    for i=1:nCells
        poly!(ax4,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=5) #:bwr
    end
    Colorbar(grid[1,2],limits=curlLims,colormap=:bwr,flipaxis=false)

    show==1 ? display(fig) : nothing
    save("$dataDirectory/pdf/internalVertexCurls.pdf",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/internalVertexCurls.pdf",fig)
    save("$dataDirectory/svg/internalVertexCurls.svg",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/internalVertexCurls.svg",fig)
    save("$dataDirectory/png/internalVertexCurls.png",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/internalVertexCurls.png",fig)
end
