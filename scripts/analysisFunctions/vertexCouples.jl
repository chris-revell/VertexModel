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

function vertexCouples(dataDirectory, show)
    isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png")
    isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf")
    isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg")
    isdir("$dataDirectory/png") ? nothing : mkpath("$dataDirectory/png")
    isdir("$dataDirectory/pdf") ? nothing : mkpath("$dataDirectory/pdf")
    isdir("$dataDirectory/svg") ? nothing : mkpath("$dataDirectory/svg")

    # Import system data
    conditionsDict    = load("$dataDirectory/dataFinal.jld2")
    @unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,L₀,A₀,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
    matricesDict = load("$dataDirectory/matricesFinal.jld2")
    @unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths,cellPressures = matricesDict["matrices"]

    T = makeCellLinks(conditionsDict["params"],matricesDict["matrices"])

    linkTriangles = makeLinkTriangles(conditionsDict["params"],matricesDict["matrices"])
    linkTriangleAreas = abs.(area.(linkTriangles))

    cellPolygons = makeCellPolygons(conditionsDict["params"],matricesDict["matrices"])

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

    lims = (-maximum(abs.(vertexCouples)),maximum(abs.(vertexCouples)))

    # Set up figure canvas
    fig = Figure(resolution=(1000,1000))
    ax = Axis(fig[1,1],aspect=DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    for k=1:nVerts
        poly!(ax,linkTriangles[k],color=[vertexCouples[k]],colorrange=lims,colormap=:bwr,strokewidth=2,strokecolor=(:black,0.0)) #:bwr
    end
    for i=1:nCells
        poly!(ax,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=5) #:bwr
    end
    Colorbar(fig[1,2],limits=lims,colormap=:bwr,flipaxis=false) #:bwr

    show==1 ? display(fig) : nothing
    save("$dataDirectory/pdf/vertexCouples.pdf",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/vertexCouples.pdf",fig)
    save("$dataDirectory/svg/vertexCouples.svg",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/vertexCouples.svg",fig)
    save("$dataDirectory/png/vertexCouples.png",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/vertexCouples.png",fig)
end
