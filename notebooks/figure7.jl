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

dataDirs = ["data/figure7/$x" for x in readdir("data/figure7") if (isdir("data/figure7/$x") && x!="old")]

labels = [L"(a)",L"(b)",L"(c)",L"(d)",L"(e)",L"(f)",L"(g)",L"(h)",L"(i)",L"(j)",L"(k)",L"(l)",L"(m)",L"(n)",L"(o)",L"(p)",L"(q)",L"(r)",L"(s)",L"(t)",L"(u)",L"(v)",L"(w)",L"(x)",L"(y)",L"(z)"]

for (i,dataDirectory) in enumerate(dataDirs)

    # Import system data
    conditionsDict    = load("$dataDirectory/dataFinal.jld2")
    @unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,L₀,A₀,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
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
    eigenvectorsLf = (eigen(Matrix(Lf))).vectors
    eigenvaluesLf = (eigen(Matrix(Lf))).values
    ḡ = ((onesVec'*H*cellDivs)/(onesVec'*H*ones(nCells))).*onesVec
    ğ = cellDivs.-ḡ
    ψ̆ = zeros(nCells)
    eigenmodeAmplitudesLf = Float64[]
    for k=2:nCells
        numerator = eigenvectorsLf[:,k]'*H*ğ
        denominator = eigenvaluesLf[k]*(eigenvectorsLf[:,k]'*H*eigenvectorsLf[:,k])
        ψ̆ .-= (numerator/denominator).*eigenvectorsLf[:,k]
        push!(eigenmodeAmplitudesLf,(numerator/denominator))
    end
    ψ̆Lims = (-maximum(abs.(ψ̆)),maximum(abs.(ψ̆)))
    ax1 = Axis(fig[i,1],aspect=DataAspect(),fontsize=32)
    hidedecorations!(ax1)
    hidespines!(ax1)
    for i=1:nCells
        poly!(ax1,cellPolygons[i],color=[ψ̆[i]],colormap=:bwr,colorrange=ψ̆Lims, strokecolor=(:black,1.0),strokewidth=1)
    end
    Label(fig[i,1,Bottom()],labels[3*i-2],textsize = 32)
    # Colorbar(fig[1,1][1,2],limits=ψ̆Lims,colormap=:bwr,flipaxis=false,align=:left)

    # psi_v axis
    Lₜ = makeLt(conditionsDict["params"],matricesDict["matrices"],T,linkTriangleAreas,trapeziumAreas)
    eigenvectorsLt = (eigen(Matrix(Lₜ))).vectors
    eigenvaluesLt = (eigen(Matrix(Lₜ))).values
    vertexDivs = -1.0.*calculateVertexDivs(conditionsDict["params"],matricesDict["matrices"],q,linkTriangleAreas)
    onesVec = ones(nVerts)
    E = Diagonal(linkTriangleAreas)
    ḡ = ((onesVec'*E*vertexDivs)/(onesVec'*E*ones(nVerts))).*onesVec
    ğ = vertexDivs.-ḡ
    ψ̆ = zeros(nVerts)
    eigenmodeAmplitudesLt = Float64[]
    for k=2:nVerts
        numerator = -eigenvectorsLt[:,k]'*E*ğ
        denominator = eigenvaluesLt[k]*(eigenvectorsLt[:,k]'*E*eigenvectorsLt[:,k])
        ψ̆ .+= (numerator/denominator).*eigenvectorsLt[:,k]
        push!(eigenmodeAmplitudesLt,(numerator/denominator))
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
    Label(fig[i,2,Bottom()],labels[3*i-1],textsize = 32)

    ax3 = Axis(fig[i,3],yscale=log10)#, xlabel="Eigenmode number", ylabel="Amplitude",fontsize=32)
    hidedecorations!(ax3)
    xlims!(ax3,0,nVerts)
    ylims!(ax3,0,1.1*maximum(abs.([eigenmodeAmplitudesLt; eigenmodeAmplitudesLf])))
    barplot!(ax3,collect(2:nCells),abs.(eigenmodeAmplitudesLf),width=1.0,color=(:blue,0.75),strokecolor=(:blue,0.75))
    barplot!(ax3,collect(2:nVerts),abs.(eigenmodeAmplitudesLt),width=1.0,color=(:orange,0.75),strokecolor=(:orange,0.75))
    Label(fig[i,3,Bottom()],labels[3*i],textsize = 32)
end

colsize!(fig.layout, 1, Aspect(1, 1))
colsize!(fig.layout, 2, Aspect(1, 1))
colsize!(fig.layout, 3, Aspect(1, 1.5))

colgap!(fig.layout,Relative(0.0))
rowgap!(fig.layout,Relative(0.01))

resize_to_layout!(fig)

display(fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/figure7.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/figure7.eps",fig)
