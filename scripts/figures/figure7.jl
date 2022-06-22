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
using UnicodeFun

# Local modules
include("$(projectdir())/scripts/analysisFunctions/functions.jl")

set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
fig = Figure(resolution=(2000,3000),fontsize = 32)

dataDirs = ["data/figure7/$x" for x in readdir("data/figure7") if (isdir("data/figure7/$x") && x!="unfinished")]

labels = [L"(a)",L"(b)",L"(c)",L"(d)",L"(e)",L"(f)",L"(g)",L"(h)",L"(i)",L"(j)",L"(k)",L"(l)",L"(m)",L"(n)",L"(o)",L"(p)",L"(q)",L"(r)",L"(s)",L"(t)",L"(u)",L"(v)",L"(w)",L"(x)",L"(y)",L"(z)"]

for (i,dataDirectory) in enumerate(dataDirs[1:5])

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
    ax1 = Axis(fig[i,1],aspect=DataAspect())
    hidedecorations!(ax1)
    hidespines!(ax1)
    for i=1:nCells
        poly!(ax1,cellPolygons[i],color=[ψ̆[i]],colormap=:bwr,colorrange=ψ̆Lims, strokecolor=(:black,1.0),strokewidth=1)
    end
    Label(fig[i,1,Bottom()],labels[3*i-2],textsize = 64)
    # Colorbar(fig[1,1][1,2],limits=ψ̆Lims,colormap=:bwr,flipaxis=false)

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
    ax2 = Axis(fig[i,2],aspect=DataAspect())
    hidedecorations!(ax2)
    hidespines!(ax2)
    for k=1:nVerts
        poly!(ax2,linkTriangles[k],color=[ψ̆[k]],colorrange=ψ̆Lims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
    end
    for i=1:nCells
        poly!(ax2,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
    end
    Label(fig[i,2,Bottom()],labels[3*i-1],textsize = 64)

    #ax3 = Axis(fig[i,3],aspect=AxisAspect(2.0), xlabel="Eigenmode number", ylabel=L"log_{10}(Amplitude)",fontsize=48,ygridvisible=false,xgridvisible=false)
    ax3 = Axis(fig[i,3],aspect=AxisAspect(2.0), xlabel="Eigenmode number", ylabel="log$(to_subscript("10"))(Amplitude)",ygridvisible=false,xgridvisible=false)

    # hidedecorations!(ax3)
    xlims!(ax3,1,nVerts)
    ylims!(ax3,-6,1)
    lines!(ax3,collect(2:nCells),log10.(abs.(eigenmodeAmplitudesLf)),linewidth=2,color=(:blue,0.5),strokecolor=(:blue,0.5))
    lines!(ax3,collect(2:nVerts),log10.(abs.(eigenmodeAmplitudesLt)),linewidth=2,color=(:orange,0.5),strokecolor=(:orange,0.5))
    Label(fig[i,3,Bottom()],labels[3*i],textsize = 64)
end

colsize!(fig.layout, 1, Aspect(1, 1))
colsize!(fig.layout, 2, Aspect(1, 1))
colsize!(fig.layout, 3, Aspect(1, 1.5))

colgap!(fig.layout,Relative(0.0))
rowgap!(fig.layout,Relative(0.01))

resize_to_layout!(fig)

display(fig)
save("$(datadir())/figure7/figure7.png",fig)
save("$(datadir())/figure7/figure7.eps",fig)