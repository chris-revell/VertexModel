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

function psicSpectrum(dataDirectory, show)
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
        ψ̆ .+= (numerator/denominator).*eigenvectors[:,k]
        push!(eigenmodeAmplitudes,(numerator/denominator))
    end

    fig = Figure(resolution=(1000,700),fontsize = 24)
    ax2 = Axis(fig[1,1], xlabel="Eigenmode number, i", ylabel="Amplitude",fontsize=32)
    xlims!(ax2,0,nCells)
    ylims!(ax2,0,1.1*maximum(abs.(eigenmodeAmplitudes)))
    barplot!(ax2,collect(2:nCells),abs.(eigenmodeAmplitudes),width=1.0,color=:blue,strokecolor=:blue)

    show==1 ? display(fig) : nothing
    save("$dataDirectory/pdf/psicSpectrum.pdf",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/psicSpectrum.pdf",fig)
    save("$dataDirectory/svg/psicSpectrum.svg",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/psicSpectrum.svg",fig)
    save("$dataDirectory/png/psicSpectrum.png",fig)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/psicSpectrum.png",fig)
end