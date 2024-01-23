# Script to produce a movie of Airy Stress evolution over vertex network for a given system

using JLD2
using SparseArrays
using LinearAlgebra
using DrWatson
using DataFrames
using FromFile
using UnPack
using CairoMakie
using Printf
using Colors

@from "$(projectdir())/src/VertexModelContainers.jl" using VertexModelContainers
@from "$(projectdir())/src/OrderAroundCell.jl" using OrderAroundCell
@from "$(projectdir())/src/AnalysisFunctions.jl" using AnalysisFunctions
@from "$(projectdir())/src/Laplacians.jl" using Laplacians
@from "$(projectdir())/src/Potentials.jl" using Potentials

for f in [f for f in readdir(datadir("sims/examples")) if occursin("γ",f)]

    folderName = "sims/examples/$f"

    fig = CairoMakie.Figure(size=(1000,1000))
    ax = Axis(fig[1,1][1,1],aspect=DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    mov = VideoStream(fig, framerate=5)

    # Skip the starting configuration (t=0)

    potentials = Vector{Float64}[]
    ψ̆Lims = [0.0,0.0]
    for t=1:99
        @unpack R, matrices, params = load(datadir(folderName,"frames","systemData$(@sprintf("%03d", t)).jld2"))
        ψ̆, spectrum = psivPotential(R,params,matrices)
        push!(potentials,ψ̆)
        ψ̆Lims .= [-max(maximum(abs.(ψ̆)),maximum(ψ̆Lims)),max(maximum(abs.(ψ̆)),maximum(ψ̆Lims))]
    end 

    Colorbar(fig[1,1][1,2],limits=ψ̆Lims,colormap=:bwr,flipaxis=true)

    for t=1:99
        @unpack R, matrices, params = load(datadir(folderName,"frames","systemData$(@sprintf("%03d", t)).jld2"))
        @unpack B, Bᵀ, C, cellPositions = matrices
        @unpack nCells, nVerts = params
        
        cellPolygons = makeCellPolygons(R,params,matrices)
        linkTriangles = makeLinkTriangles(R,params,matrices)
        empty!(ax)
        # ax.title = "t = $(@sprintf("%.2f", t))"
        for k=1:nVerts
            poly!(ax,linkTriangles[k],color=[potentials[t][k]],colorrange=ψ̆Lims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25))
        end
        for i=1:nCells
            poly!(ax,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=2)
        end
        reset_limits!(ax)
        recordframe!(mov)
    end

    save(datadir(folderName,"movieAiryStressV.mp4"),mov)

end