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

function psicPotential(dataDirectory, show)
    
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

    Lf = makeLf(conditionsDict["params"],matricesDict["matrices"],trapeziumAreas)

    # Calculate div on each cell
    cellDivs = -1.0.*calculateCellDivs(conditionsDict["params"],matricesDict["matrices"])

    onesVec = ones(nCells)
    H = Diagonal(cellAreas)

    eigenvectors = (eigen(Matrix(Lf))).vectors
    eigenvalues = (eigen(Matrix(Lf))).values

    ḡ = ((onesVec'*H*cellDivs)/(onesVec'*H*ones(nCells))).*onesVec
    ğ = cellDivs.-ḡ
    ψ̆ = zeros(nCells)
    eigenmodeAmplitudes = Float64[]
    for k=2:nCells
        numerator = eigenvectors[:,k]'*H*ğ
        denominator = eigenvalues[k]*(eigenvectors[:,k]'*H*eigenvectors[:,k])
        ψ̆ .-= (numerator/denominator).*eigenvectors[:,k]
        push!(eigenmodeAmplitudes,(numerator/denominator))
    end

    ψ̆Lims = (-maximum(abs.(ψ̆)),maximum(abs.(ψ̆)))

    fig = Figure(resolution=(1000,1000),fontsize = 32)
    ax1 = Axis(fig[1,1][1,1],aspect=DataAspect())
    hidedecorations!(ax1)
    hidespines!(ax1)
    for i=1:nCells
        poly!(ax1,cellPolygons[i],color=[ψ̆[i]],colormap=:bwr,colorrange=ψ̆Lims, strokecolor=(:black,1.0),strokewidth=5)
    end
    Colorbar(fig[1,1][1,2],limits=ψ̆Lims,colormap=:bwr,flipaxis=false)

    show==1 ? display(fig) : nothing
    # save("$dataDirectory/pdf/psicPotential.pdf",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/psicPotential.pdf",fig)
    # save("$dataDirectory/svg/psicPotential.svg",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/psicPotential.svg",fig)
    save("$dataDirectory/psicPotential.png",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/psicPotential.png",fig)
end

psicPotential(datadir("annealing","L₀=3.0_γ=0.15_23-02-08-08-03-02"),1)