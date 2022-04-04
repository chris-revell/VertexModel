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

function psivSpectrum(dataDirectory, show)
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


    vertexDivs = -1.0.*calculateVertexDivs(conditionsDict["params"],matricesDict["matrices"],T,linkTriangleAreas)

    onesVec = ones(nVerts)
    E = Diagonal(linkTriangleAreas)

    ḡ = ((onesVec'*E*vertexDivs)/(onesVec'*E*ones(nVerts))).*onesVec
    ğ = vertexDivs.-ḡ
    ψ̆ = zeros(nVerts)
    eigenmodeAmplitudes = Float64[]
    for k=2:nVerts
        numerator = -eigenvectors[:,k]'*E*ğ
        denominator = eigenvalues[k]*(eigenvectors[:,k]'*E*eigenvectors[:,k])
        ψ̆ .+= (numerator/denominator).*eigenvectors[:,k]
        push!(eigenmodeAmplitudes,(numerator/denominator))
    end

    fig = Figure(resolution=(1000,700),fontsize = 24)
    ax2 = Axis(fig[1,1], xlabel="Eigenmode number, k", ylabel="Amplitude",fontsize=32)
    xlims!(ax2,0,nVerts)
    ylims!(ax2,0,1.1*maximum(abs.(eigenmodeAmplitudes)))
    barplot!(ax2,collect(2:nVerts),abs.(eigenmodeAmplitudes),width=1.0,color=:blue,strokecolor=:blue)

    show==1 ? display(fig) : nothing
    save("$dataDirectory/pdf/psivSpectrum.pdf",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/psivSpectrum.pdf",fig)
    save("$dataDirectory/svg/psivSpectrum.svg",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/psivSpectrum.svg",fig)
    save("$dataDirectory/png/psivSpectrum.png",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/psivSpectrum.png",fig)
end
