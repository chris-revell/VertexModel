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
includet("$(projectdir())/notebooks/functions.jl")

function laplacianTableauLf(dataDirectory, show)

    isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png")
    isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf")
    isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg")

    # Import system data
    conditionsDict    = load("$dataDirectory/dataFinal.jld2")
    @unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
    matricesDict = load("$dataDirectory/matricesFinal.jld2")
    @unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths = matricesDict["matrices"]

    T = makeCellLinks(conditionsDict["params"],matricesDict["matrices"])

    edgeTrapezia = makeEdgeTrapezia(conditionsDict["params"],matricesDict["matrices"])
    trapeziumAreas = abs.(area.(edgeTrapezia))

    linkTriangles = makeLinkTriangles(conditionsDict["params"],matricesDict["matrices"])
    linkTriangleAreas = abs.(area.(linkTriangles))

    cellPolygons = makeCellPolygons(conditionsDict["params"],matricesDict["matrices"])

    Lf = makeLf(conditionsDict["params"],matricesDict["matrices"],trapeziumAreas)

    decomposition = (eigen(Matrix(Lf))).vectors

    signInversions = [3,4,5,9,11,12,16,17,18,19,21]

    # Set up figure canvas
    fig = Figure(resolution=(750,1600))
    grid = fig[1,1] = GridLayout()
    for x=1:5
        for y=1:4
            eigenvectorIndex = ((y-1)*5 + x)+1
            lims = (-maximum(abs.(decomposition[:,eigenvectorIndex])),maximum(abs.(decomposition[:,eigenvectorIndex])))
            ax = Axis(grid[y,x],aspect=DataAspect())
            hidedecorations!(ax)
            hidespines!(ax)
            if eigenvectorIndex in signInversions
                for i=1:nCells
                    poly!(ax,cellPolygons[i],color=[-decomposition[i,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokecolor=(:black,1.0),strokewidth=1) #:bwr
                end
            else
                for i=1:nCells
                    poly!(ax,cellPolygons[i],color=[decomposition[i,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokecolor=(:black,1.0),strokewidth=1) #:bwr
                end
            end
            Label(grid[y,x,Bottom()],
                    L"i=%$eigenvectorIndex",
                    textsize = 16,
            )
        end
    end
    for x=1:5
        for y=1:4
            eigenvectorIndex = ((y-1)*5 + x)+(nCells-20)
            lims = (-maximum(abs.(decomposition[:,eigenvectorIndex])),maximum(abs.(decomposition[:,eigenvectorIndex])))
            ax = Axis(grid[y+4,x],aspect=DataAspect())
            hidedecorations!(ax)
            hidespines!(ax)
            for i=1:nCells
                poly!(ax,cellPolygons[i],color=[decomposition[i,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokecolor=(:black,1.0),strokewidth=1) #:bwr
            end
            Label(grid[y+4,x,Bottom()],
                    L"i=%$eigenvectorIndex",
                    textsize = 16,
            )
        end
    end


    show==1 ? display(fig) : nothing
    save("$dataDirectory/eigenvectorTableauLf.pdf",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/eigenvectorTableauLf.pdf",fig)
    save("$dataDirectory/eigenvectorTableauLf.svg",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/eigenvectorTableauLf.svg",fig)
    save("$dataDirectory/eigenvectorTableauLf.png",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/eigenvectorTableauLf.png",fig)
end
