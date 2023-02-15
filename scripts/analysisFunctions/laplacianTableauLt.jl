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
using FromFile

# Local modules
@from "$(projectdir())/src/AnalysisFunctions.jl" using AnalysisFunctions 

function laplacianTableauLt(dataDirectory, show)

    # isdir("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/png") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/png")
    # isdir("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf")
    # isdir("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg")
    # isdir("$dataDirectory/png") ? nothing : mkpath("$dataDirectory/png")
    # isdir("$dataDirectory/pdf") ? nothing : mkpath("$dataDirectory/pdf")
    # isdir("$dataDirectory/svg") ? nothing : mkpath("$dataDirectory/svg")

    # Import system data
    conditionsDict    = load("$dataDirectory/dataFinal.jld2")
    @unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
    matricesDict = load("$dataDirectory/matricesFinal.jld2")
    @unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths = matricesDict["matrices"]

    T = makeCellLinks(conditionsDict["params"],matricesDict["matrices"])

    edgeTrapezia = makeEdgeTrapezia(conditionsDict["params"],matricesDict["matrices"])
    trapeziumAreas = abs.(area.(edgeTrapezia))

    linkTriangles = makeLinkTriangles(conditionsDict["params"],matricesDict["matrices"])
    linkTriangleAreas = abs.(area.(linkTriangles))

    cellPolygons = makeCellPolygons(conditionsDict["params"],matricesDict["matrices"])

    Lₜ = makeLt(conditionsDict["params"],matricesDict["matrices"],T,linkTriangleAreas,trapeziumAreas)

    decomposition = (eigen(Matrix(Lₜ))).vectors

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
            for k=1:nVerts
                poly!(ax,linkTriangles[k],color=[decomposition[k,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
            end
            for i=1:nCells
                poly!(ax,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
            end
            Label(grid[y,x,Bottom()],
                    L"k=%$eigenvectorIndex",
                    fontsize = 16,
            )
        end
    end
    for x=1:5
        for y=1:4
            eigenvectorIndex = ((y-1)*5 + x)+(nVerts-20)
            lims = (-maximum(abs.(decomposition[:,eigenvectorIndex])),maximum(abs.(decomposition[:,eigenvectorIndex])))
            ax = Axis(grid[y+4,x],aspect=DataAspect())
            hidedecorations!(ax)
            hidespines!(ax)
            for k=1:nVerts
                poly!(ax,linkTriangles[k],color=[decomposition[k,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
            end
            for i=1:nCells
                poly!(ax,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
            end
            Label(grid[y+4,x,Bottom()],
                    L"k=%$eigenvectorIndex",
                    fontsize = 16,
            )
        end
    end

    show==1 ? display(fig) : nothing
    # save("$dataDirectory/svg/eigenvectorTableauLt.svg",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/eigenvectorTableauLt.svg",fig)
    # save("$dataDirectory/pdf/eigenvectorTableauLt.pdf",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/eigenvectorTableauLt.pdf",fig)
    save("$dataDirectory/eigenvectorTableauLt.png",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/eigenvectorTableauLt.png",fig)
end

dataDirectory = datadir("annealing","L₀=0.5_γ=0.1_23-02-02-18-06-16")
laplacianTableauLt(dataDirectory, 1)