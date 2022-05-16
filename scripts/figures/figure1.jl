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
includet("$(projectdir())/scripts/analysisFunctions/functions.jl")

# dataDirectory = "data/old/2022-02-28-19-30-22"
#
# centralCell = 14
#
# plotCells = 1
# plotLinks = 1
# scatterVertices = 0
# scatterEdges = 0
# scatterCells = 0
# annotateVertices = 0
# annotateEdges = 0
# annotateCells = 0

function figure1(dataDirectory,centralCell,plotCells,plotLinks,scatterVertices,scatterEdges,scatterCells,annotateVertices,annotateEdges,annotateCells)

    # Import system data
    conditionsDict    = load("$dataDirectory/dataFinal.jld2")
    @unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
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
    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    fig = Figure(resolution=(1000,400))
    ax1 = Axis(fig[1,1],aspect=DataAspect())
    hidedecorations!(ax1)
    hidespines!(ax1)
    Label(fig[1,1,Bottom()],L"(a)",textsize = 32)

    if plotCells==1
        # Plot cell polygons
        for i=1:nCells
            if cellNeighbourMatrix[centralCell,i] == 0
                poly!(ax1,cellPolygons[i],color=(getRandomColor(i),0.25),strokecolor=(:black,0.5),strokewidth=1)
            else
                poly!(ax1,cellPolygons[i],color=(getRandomColor(i),1.0),strokecolor=(:black,1.0),strokewidth=2)
            end
        end
    end

    if plotLinks==1
        for j=1:nEdges
            edgeCells = findall(!iszero,B[:,j])
            if boundaryEdges[j] == 0
                lines!(ax1,Point2f.(cellPositions[edgeCells]),linewidth=1,color=(:white,1.0))
            else
                lines!(ax1,Point2f.([cellPositions[edgeCells[1]],cᵖ[j]]),linewidth=1,color=(:white,1.0))
            end
        end
    end

    if scatterVertices==1
        # # Scatter vertex locations
        scatter!(ax1,Point2f.(R),alpha=0.5,color=:blue)
        if annotateVertices==1
            annotations!(ax1,string.(collect(1:nVerts)),Point2f.(R),color=:blue)
        end
    end

    # Edge labels
    if scatterEdges==1
        scatter!(ax1,Point2f.(edgeMidpoints),alpha=0.5,color=:green)
        if annotateEdges==1
            annotations!(ax1,string.(collect(1:nEdges)),Point2f.(edgeMidpoints),color=:green)
        end
    end

    # Scatter cell centroid locations
    if scatterCells==1
        scatter!(ax1,Point2f.(cellPositions),color=:red)
        if annotateCells==1
            annotations!(ax1,string.(collect(1:nCells)),Point2f.(cellPositions),color=:red)
        end
    end

    ax2 = Axis(fig[1,2],aspect=DataAspect())
    hidedecorations!(ax2)
    hidespines!(ax2)
    Label(fig[1,2,Bottom()],L"(b)",textsize = 32)
    image!(ax2,rotr90(load("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/OtherFigures/geometryCb.png")))

    ax3 = Axis(fig[1,3],aspect=DataAspect())
    hidedecorations!(ax3)
    hidespines!(ax3)
    Label(fig[1,3,Bottom()],L"(c)",textsize = 32)
    image!(ax3,rotr90(load("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/OtherFigures/geometryCc.png")))

    # boxes = [Box(g,color=(:white,0.0)) for g in [fig[1,1],fig[1,2],fig[1,3]]]

    colgap!(fig.layout,1,Relative(0.0))
    colgap!(fig.layout,2,Relative(0.03))
    resize_to_layout!(fig)

    # display(fig)
    save("$dataDirectory/figure1.pdf",fig)
    # #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/figure1.pdf",fig)
    # save("$dataDirectory/figure1.svg",fig)
    # #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/figure1.svg",fig)
    save("$dataDirectory/figure1.png",fig)
    # #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/figure1.eps",fig)
    # #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/figure1.png",fig)
end
