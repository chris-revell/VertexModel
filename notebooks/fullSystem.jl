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
includet("$(projectdir())/notebooks/functions.jl")

# centralCell=1 #14

function fullSystem(dataDirectory, centralCell, plotCells, plotLinks, scatterEdges, annotateEdges, scatterVertices, annotateVertices, scatterCells, annotateCells, show)
    isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png")
    isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf")
    isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg")
    isdir("$dataDirectory/png") ? nothing : mkpath("$dataDirectory/png")
    isdir("$dataDirectory/pdf") ? nothing : mkpath("$dataDirectory/pdf")
    isdir("$dataDirectory/svg") ? nothing : mkpath("$dataDirectory/svg")

    # Import system data
    conditionsDict    = load("$dataDirectory/dataFinal.jld2")
    @unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
    matricesDict = load("$dataDirectory/matricesFinal.jld2")
    @unpack A,B,Bᵀ,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ = matricesDict["matrices"]

    cellNeighbourMatrix = B*Bᵀ
    dropzeros!(cellNeighbourMatrix)

    onesVec = ones(1,nCells)
    boundaryEdges = abs.(onesVec*B)
    cᵖ = boundaryEdges'.*edgeMidpoints

    # Find cell midpoint links T
    T = makeCellLinks(conditionsDict["params"],matricesDict["matrices"])

    # Create vector of polygons for each cellwparams"],matricesDict["matrices"])
    cellPolygons = makeCellPolygons(conditionsDict["params"],matricesDict["matrices"])

    # Set up figure canvas
    fig = Figure(resolution=(1000,1000))
    ax = Axis(fig[1,1],aspect=DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)

    if plotCells==1
        # Plot cell polygons
        for i=1:nCells
            if cellNeighbourMatrix[centralCell,i] == 0
                poly!(ax,cellPolygons[i],color=(getRandomColor(i),0.25),strokecolor=(:black,0.5),strokewidth=1)
            else
                poly!(ax,cellPolygons[i],color=(getRandomColor(i),1.0),strokecolor=(:black,1.0),strokewidth=2)
            end
        end
    end

    if plotLinks==1
        for j=1:nEdges
            edgeCells = findall(!iszero,B[:,j])
            if boundaryEdges[j] == 0
                lines!(ax,Point2f.(cellPositions[edgeCells]),linewidth=1,color=(:white,1.0))
            else
                lines!(ax,Point2f.([cellPositions[edgeCells[1]],cᵖ[j]]),linewidth=1,color=(:white,1.0))
            end
        end
    end

    if scatterVertices==1
        # # Scatter vertex locations
        scatter!(ax,Point2f.(R),alpha=0.5,color=:blue)
        if annotateVertices==1
            annotations!(ax,string.(collect(1:nVerts)),Point2f.(R),color=:blue)
        end
    end

    # Edge labels
    if scatterEdges==1
        scatter!(ax,Point2f.(edgeMidpoints),alpha=0.5,color=:green)
        if annotateEdges==1
            annotations!(ax,string.(collect(1:nEdges)),Point2f.(edgeMidpoints),color=:green)
        end
    end

    # Scatter cell centroid locations
    if scatterCells==1
        scatter!(ax,Point2f.(cellPositions),color=:red)
        if annotateCells==1
            annotations!(ax,string.(collect(1:nCells)),Point2f.(cellPositions),color=:red)
        end
    end

    show==1 ? display(fig) : nothing
    save("$dataDirectory/pdf/fullSystemCentralCell$centralCell.pdf",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/fullSystemCentralCell$centralCell.pdf",fig)
    save("$dataDirectory/svg/fullSystemCentralCell$centralCell.svg",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/fullSystemCentralCell$centralCell.svg",fig)
    save("$dataDirectory/png/fullSystemCentralCell$centralCell.png",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/fullSystemCentralCell$centralCell.png",fig)
end
