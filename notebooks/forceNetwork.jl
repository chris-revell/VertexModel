# Import Julia packages
using DrWatson
@quickactivate
using Revise
using LinearAlgebra
using SparseArrays
using StaticArrays
using CairoMakie
using UnPack
using GeometryBasics
using DelimitedFiles
using Random
using Colors
using JLD2

# Local modules
includet("$(projectdir())/notebooks/functions.jl")

# centralCell=1 #14

function forceNetwork(dataDirectory, centralCell, show)
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
    @unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeLengths,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths = matricesDict["matrices"]

    # Set up figure canvas
    fig = Figure(resolution=(1000,1000))
    ax2 = Axis(fig[1,1],aspect=DataAspect())
    hidedecorations!(ax2)
    hidespines!(ax2)

    # Find all cells neighbouring original cell
    cellNeighbourMatrix = B*Bᵀ
    dropzeros!(cellNeighbourMatrix)
    neighbouringCells = findall(!iszero,cellNeighbourMatrix[centralCell,:])

    # Find and sort all vertices around cells neighbouring centralCell
    cellVerticesDict = makeCellVerticesDict(conditionsDict["params"],matricesDict["matrices"])

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
    startPosition = [SVector{2,Float64}([0.0,0.0])] #Make this a single component array of SVectors to avoid scope issues with immutable objects
    for (i,v) in enumerate(cellVerticesDict[centralCell])
        arrows!(ax2,Point2f.([startPosition[1]]),Vec2f.([ϵ*F[v,centralCell]]),linewidth=6,arrowsize=30,color=(getRandomColor(centralCell),0.9))
        startPosition[1] = startPosition[1] + ϵ*F[v,centralCell]
        H = Array{SVector{2,Float64}}(undef,length(cellVerticesDict[neighbouringCells[i]])+1)
        cellForces = SVector{2, Float64}[]
        # Circular permutation of vertices to ensure vertex v is the first index
        # in the ordered cellVertices list around cell neighbouringCells[i]
        index = findall(x->x==v, cellVerticesDict[neighbouringCells[i]])
        cellVertices = circshift(cellVerticesDict[neighbouringCells[i]],1-index[1])
        H[1] = startPosition[1]
        for (j,cv) in enumerate(cellVertices)
            push!(cellForces,+ϵ*F[cv,neighbouringCells[i]])
            H[j+1] = H[j]+cellForces[end]
        end
        arrows!(ax2,Point2f.(H),Vec2f.(cellForces),color=(getRandomColor(neighbouringCells[i]),0.9),linewidth=6,arrowsize=30)
    end

    show==1 ? display(fig) : nothing

    save("$dataDirectory/pdf/cell$(centralCell)ForceNetwork.pdf",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/cell$(centralCell)ForceNetwork.pdf",fig)
    save("$dataDirectory/svg/cell$(centralCell)ForceNetwork.svg",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/cell$(centralCell)ForceNetwork.svg",fig)
    save("$dataDirectory/png/cell$(centralCell)ForceNetwork.png",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/cell$(centralCell)ForceNetwork.png",fig)
end
