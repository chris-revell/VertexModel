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

folderName = "newlongTest/L₀=0.75_realTimetMax=86400.0_t1Threshold=0.01_γ=0.2_23-03-08-20-49-23"

fig = CairoMakie.Figure(resolution=(1000,1000))
ax = Axis(fig[1,1][1,1],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
mov = VideoStream(fig, framerate=5)

potentials = Vector{Float64}[]
ψ̆Lims = [0.0,0.0]
for t=0:100
    @unpack R, matrices, params = load(datadir(folderName,"frames","systemData$(@sprintf("%03d", t)).jld2"))
    ψ̆, spectrum = psicPotential(R,params,matrices)
    push!(potentials,ψ̆)
    ψ̆Lims .= [-max(maximum(abs.(ψ̆)),maximum(ψ̆Lims)),max(maximum(abs.(ψ̆)),maximum(ψ̆Lims))]
end 

Colorbar(fig[1,1][1,2],limits=ψ̆Lims,colormap=:bwr,flipaxis=true)


for t=0:100
    @unpack R, matrices, params = load(datadir(folderName,"frames","systemData$(@sprintf("%03d", t)).jld2"))
    @unpack B, Bᵀ, C, cellPositions = matrices
    @unpack nCells = params

    cellPolygons = makeCellPolygons(R,params,matrices)
    empty!(ax)
    ax.title = "t = $(@sprintf("%.2f", t))"
    xlims!(ax, 1.1*min(minimum(first.(R)), minimum(last.(R))), 1.1*max(maximum(first.(R)), maximum(last.(R))))
    ylims!(ax, 1.1*min(minimum(first.(R)), minimum(last.(R))), 1.1*max(maximum(first.(R)), maximum(last.(R))))
    for i=1:nCells
        poly!(ax,cellPolygons[i],color=[potentials[t+1][i]],colormap=:bwr,colorrange=Tuple(ψ̆Lims), strokecolor=(:black,1.0),strokewidth=5)
    end
    recordframe!(mov)
end

save(datadir(folderName,"airyStressC.mp4"),mov)
