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

function capitalPsivPotentialLtDerivative(dataDirectory,show)
    isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png")
    isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf")
    isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg")
    isdir("$dataDirectory/png") ? nothing : mkpath("$dataDirectory/png")
    isdir("$dataDirectory/pdf") ? nothing : mkpath("$dataDirectory/pdf")
    isdir("$dataDirectory/svg") ? nothing : mkpath("$dataDirectory/svg")

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

    eigenvectors = (eigen(Matrix(Lₜ))).vectors
    eigenvalues = (eigen(Matrix(Lₜ))).values

    q = calculateSpokes(conditionsDict["params"],matricesDict["matrices"])

    vertexCurls = calculateVertexCurls(conditionsDict["params"],matricesDict["matrices"],q,linkTriangleAreas)

    onesVec = ones(nVerts)
    E = Diagonal(linkTriangleAreas)

    ḡ = ((onesVec'*E*vertexCurls)/(onesVec'*E*ones(nVerts))).*onesVec
    ğ = vertexCurls.-ḡ
    ψ̆ = zeros(nVerts)
    eigenmodeAmplitudes = Float64[]
    for k=2:nVerts
        numerator = -eigenvectors[:,k]'*E*ğ
        denominator = eigenvalues[k]*(eigenvectors[:,k]'*E*eigenvectors[:,k])
        ψ̆ .+= (numerator/denominator).*eigenvectors[:,k]
        push!(eigenmodeAmplitudes,(numerator/denominator))
    end

    ψ̆Lims = (-maximum(abs.(ψ̆)),maximum(abs.(ψ̆)))

    derivative = Lₜ*ψ̆

    derivativeLims = (-maximum(abs.(derivative)),maximum(abs.(derivative)))

    fig = Figure(resolution=(1000,1000),fontsize = 24)
    ax1 = Axis(fig[1,1][1,1],aspect=DataAspect(),fontsize=32)
    hidedecorations!(ax1)
    hidespines!(ax1)
    for k=1:nVerts
        poly!(ax1,linkTriangles[k],color=[derivative[k]],colorrange=derivativeLims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
    end
    for i=1:nCells
        poly!(ax1,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
    end
    Colorbar(fig[1,1][1,2],limits=derivativeLims,colormap=:bwr,flipaxis=false,align=:left)


    show==1 ? display(fig) : nothing
    save("$dataDirectory/pdf/capitalPsivPotentialLtDerivative.pdf",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/capitalPsivPotentialLtDerivative.pdf",fig)
    save("$dataDirectory/svg/capitalPsivPotentialLtDerivative.svg",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/capitalPsivPotentialLtDerivative.svg",fig)
    save("$dataDirectory/png/capitalPsivPotentialLtDerivative.png",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/capitalPsivPotentialLtDerivative.png",fig)
end
capitalPsivPotentialLtDerivative
