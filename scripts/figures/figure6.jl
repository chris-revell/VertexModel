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

# dataDirectory = "data/old/2022-02-28-19-30-22"

function figure6(dataDirectory)
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

    # Set up figure canvas
    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica",fontsize=19)
    fig = Figure(resolution=(1000,1400))

    # psi_c potential axis
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
    ax1 = Axis(fig[1,1],aspect=DataAspect())
    hidedecorations!(ax1)
    hidespines!(ax1)
    for i=1:nCells
        poly!(ax1,cellPolygons[i],color=[ψ̆[i]],colormap=:bwr,colorrange=ψ̆Lims, strokecolor=(:black,1.0),strokewidth=1)
    end
    Label(fig[1,1,Bottom()],L"(a)",textsize = 32)

    # psi_c spectrum axis
    ax2 = Axis(fig[1,2], xlabel="Eigenmode number, i", ylabel="Amplitude", alignmode = Outside())
    xlims!(ax2,0,nCells)
    ylims!(ax2,0,1.1*maximum(abs.(eigenmodeAmplitudes)))
    barplot!(ax2,collect(2:nCells),abs.(eigenmodeAmplitudes),width=1.0,color=:blue,strokecolor=:blue)
    Label(fig[1,2,Bottom()],L"(b)",textsize = 32)

    # psi_v potential axis
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
    ax3 = Axis(fig[2,1],aspect=DataAspect())
    hidedecorations!(ax3)
    hidespines!(ax3)
    for k=1:nVerts
        poly!(ax3,linkTriangles[k],color=[ψ̆[k]],colorrange=ψ̆Lims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
    end
    for i=1:nCells
        poly!(ax3,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
    end
    # Colorbar(fig[1,1][1,2],limits=ψ̆Lims,colormap=:bwr,flipaxis=false)
    Label(fig[2,1,Bottom()],L"(c)",textsize = 32)

    ax4 = Axis(fig[2,2], xlabel="Eigenmode number, k", ylabel="Amplitude", alignmode = Outside())
    xlims!(ax4,0,nVerts)
    ylims!(ax4,0,1.1*maximum(abs.(eigenmodeAmplitudes)))
    barplot!(ax4,collect(2:nVerts),abs.(eigenmodeAmplitudes),width=1.0,color=:orange,strokecolor=:orange)
    Label(fig[2,2,Bottom()],L"(d)",textsize = 32)

    # Capital psi_v potential axis
    Lₜ = makeLt(conditionsDict["params"],matricesDict["matrices"],T,linkTriangleAreas,trapeziumAreas)
    eigenvectors = (eigen(Matrix(Lₜ))).vectors
    eigenvalues = (eigen(Matrix(Lₜ))).values
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
    ax5 = Axis(fig[3,1],aspect=DataAspect())
    hidedecorations!(ax5)
    hidespines!(ax5)
    for k=1:nVerts
        poly!(ax5,linkTriangles[k],color=[ψ̆[k]],colorrange=ψ̆Lims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
    end
    for i=1:nCells
        poly!(ax5,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
    end
    # Colorbar(fig[1,1][1,2],limits=ψ̆Lims,colormap=:bwr,flipaxis=false)
    Label(fig[3,1,Bottom()],L"(e)",textsize = 32)

    # Capital psi_v spectrum axis
    ax6 = Axis(fig[3,2], xlabel="Eigenmode number, k", ylabel="Amplitude", alignmode = Outside())
    xlims!(ax6,0,nVerts)
    ylims!(ax6,0,1.1*maximum(abs.(eigenmodeAmplitudes)))
    barplot!(ax6,collect(2:nVerts),abs.(eigenmodeAmplitudes),width=1.0,color=:green,strokecolor=:green)
    Label(fig[3,2,Bottom()],L"(f)",textsize = 32)



    # for i=1:3
    #     for j = 1:2
    #         Box(fig[i,j],color=(:white,0.0))
    #     end
    # end

    # force the layout cell [1, 1] to be square
    colsize!(fig.layout, 1, Aspect(1, 1))
    colsize!(fig.layout, 2, Aspect(1,1.5))

    colgap!(fig.layout,Relative(0.0))
    rowgap!(fig.layout,Relative(0.01))

    resize_to_layout!(fig)

    # display(fig)
    save("$dataDirectory/figure6.png",fig)
    save("$dataDirectory/figure6.eps",fig)
end
