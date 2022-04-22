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
include("$(projectdir())/scripts/analysisFunctions/functions.jl")

fig = Figure(resolution=(500,2000),fontsize = 24)

dataDirs = ["data/sims/$x" for x in readdir("data/sims/") if isdir("data/sims/$x")]

labels = [L"(a)",L"(b)",L"(c)",L"(d)",L"(e)",L"(f)",L"(g)",L"(h)",L"(i)",L"(j)",L"(k)",L"(l)",L"(m)",L"(n)",L"(o)",L"(p)",L"(q)",L"(r)",L"(s)",L"(t)",L"(u)",L"(v)",L"(w)",L"(x)",L"(y)",L"(z)"]

for (i,dataDirectory) in enumerate(dataDirs[end-6:end])

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

    q = calculateSpokes(conditionsDict["params"],matricesDict["matrices"])

    # psi_c axis
    Lf = makeLf(conditionsDict["params"],matricesDict["matrices"],trapeziumAreas)
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
    ax1 = Axis(fig[i,1],aspect=DataAspect(),fontsize=32)
    hidedecorations!(ax1)
    hidespines!(ax1)
    for i=1:nCells
        poly!(ax1,cellPolygons[i],color=[ψ̆[i]],colormap=:bwr,colorrange=ψ̆Lims, strokecolor=(:black,1.0),strokewidth=1)
    end
    Label(fig[i,1,Bottom()],labels[2*i-1],textsize = 32)
    # Colorbar(fig[1,1][1,2],limits=ψ̆Lims,colormap=:bwr,flipaxis=false,align=:left)

    # psi_v axis
    Lₜ = makeLt(conditionsDict["params"],matricesDict["matrices"],T,linkTriangleAreas,trapeziumAreas)
    eigenvectors = (eigen(Matrix(Lₜ))).vectors
    eigenvalues = (eigen(Matrix(Lₜ))).values
    vertexDivs = -1.0.*calculateVertexDivs(conditionsDict["params"],matricesDict["matrices"],q,linkTriangleAreas)
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
    ψ̆Lims = (-maximum(abs.(ψ̆)),maximum(abs.(ψ̆)))
    ax2 = Axis(fig[i,2],aspect=DataAspect(),fontsize=32)
    hidedecorations!(ax2)
    hidespines!(ax2)
    for k=1:nVerts
        poly!(ax2,linkTriangles[k],color=[ψ̆[k]],colorrange=ψ̆Lims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
    end
    for i=1:nCells
        poly!(ax2,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
    end
    Label(fig[i,2,Bottom()],labels[2*i],textsize = 32)
end

colsize!(fig.layout, 1, Aspect(1, 1))
colsize!(fig.layout, 2, Aspect(1, 1))

colgap!(fig.layout,Relative(0.0))
rowgap!(fig.layout,Relative(0.0))

resize_to_layout!(fig)

display(fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/figure7.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/figure7.eps",fig)
