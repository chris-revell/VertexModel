# Script to produce a movie of Airy Stress evolution over cell network for a given system

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

# for f in [f for f in readdir(datadir("sims/examples")) if occursin("γ",f)]

    # folderName = "sims/examples/$f"
    folderName = "sims/L₀=0.75_nCells=800_realTimetMax=173000.0_γ=0.2_24-01-31-14-40-29"

    fig = CairoMakie.Figure(size=(1000,1000))
    ax = Axis(fig[1,1][1,1],aspect=DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    mov = VideoStream(fig, framerate=5)

    # Skip the starting configuration (t=0)

    potentials = Vector{Float64}[]
    ψ̆Lims = [0.0,0.0]
    for t=0:100
        @unpack R, matrices, params = load(datadir(folderName,"frameData","systemData$(@sprintf("%03d", t)).jld2"))
        ψ̆, spectrum = psicPotential(R,params,matrices)
        push!(potentials,ψ̆)
        ψ̆Lims .= [-max(maximum(abs.(ψ̆)),maximum(ψ̆Lims)),max(maximum(abs.(ψ̆)),maximum(ψ̆Lims))]
    end 

    Colorbar(fig[1,1][1,2],limits=ψ̆Lims,colormap=:bwr,flipaxis=true)

    for t=0:100
        @unpack R, matrices, params = load(datadir(folderName,"frameData","systemData$(@sprintf("%03d", t)).jld2"))
        @unpack B, Bᵀ, C, cellPositions = matrices
        @unpack nCells = params

        cellPolygons = makeCellPolygons(R,params,matrices)
        empty!(ax)
        # ax.title = "t = $(@sprintf("%.2f", t))"
        for i=1:nCells
            poly!(ax,cellPolygons[i],color=potentials[t+1][i],colormap=:bwr,colorrange=Tuple(ψ̆Lims), strokecolor=(:black,1.0),strokewidth=2)
        end
        reset_limits!(ax)
        recordframe!(mov)
    end

    save(datadir(folderName,"movieAiryStressC.mp4"),mov)

# end