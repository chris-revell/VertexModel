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

function psicPotentialLfDerivative(dataDirectory, show)
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


    derivative = Lf*ψ̆

    derivativeLims = (-maximum(abs.(derivative)),maximum(abs.(derivative)))

    fig = Figure(resolution=(1000,1000),fontsize = 32)
    ax1 = Axis(fig[1,1][1,1],aspect=DataAspect())
    # ax1.title = LaTeXString("L_f\phi\breve")
    hidedecorations!(ax1)
    hidespines!(ax1)
    for i=1:nCells
        poly!(ax1,cellPolygons[i],color=[derivative[i]],colormap=:bwr,colorrange=derivativeLims, strokecolor=(:black,1.0),strokewidth=5)
    end
    Colorbar(fig[1,1][1,2],limits=derivativeLims,colormap=:bwr,flipaxis=false)

    ğLims = (-maximum(abs.(ğ)),maximum(abs.(ğ)))
    fig2 = Figure(resolution=(1000,1000),fontsize = 32)
    ax2 = Axis(fig2[1,1][1,1],aspect=DataAspect())
    # ax1.title = LaTeXString("L_f\phi\breve")
    hidedecorations!(ax2)
    hidespines!(ax2)
    for i=1:nCells
        poly!(ax2,cellPolygons[i],color=[ğ[i]],colormap=:bwr,colorrange=ğLims, strokecolor=(:black,1.0),strokewidth=5)
    end
    Colorbar(fig2[1,1][1,2],limits=ğLims,colormap=:bwr,flipaxis=false)

    show==1 ? display(fig) : nothing
    save("$dataDirectory/pdf/psicPotentialLfDerivative.pdf",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/psicPotentialLfDerivative.pdf",fig)
    save("$dataDirectory/svg/psicPotentialLfDerivative.svg",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/psicPotentialLfDerivative.svg",fig)
    save("$dataDirectory/png/psicPotentialLfDerivative.png",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/psicPotentialLfDerivative.png",fig)

    # save("$dataDirectory/pdf/gbreve.pdf",fig2)
    # #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/gbreve.pdf",fig2)
    # save("$dataDirectory/svg/gbreve.svg",fig2)
    # #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/gbreve.svg",fig2)
    # save("$dataDirectory/png/gbreve.png",fig2)
    # #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/gbreve.png",fig2)
end
