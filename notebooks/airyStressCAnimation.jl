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

folderName = "annealing/L₀=3.0_γ=0.15_23-02-02-17-13-02"

fig = CairoMakie.Figure(resolution=(1000,1000))
ax = Axis(fig[1,1][1,1],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
mov = VideoStream(fig, framerate=5)

potentials = Vector{Float64}[]
ψ̆Lims = [0.0,0.0]
for t=0:100
    @unpack matrices = load(datadir(folderName,"frames","matrices$(@sprintf("%03d", t)).jld2"))
    @unpack params = load(datadir(folderName,"frames","params$(@sprintf("%03d", t)).jld2"))
    ψ̆, spectrum = psicPotential(params,matrices)
    push!(potentials,ψ̆)
    ψ̆Lims .= [-max(maximum(abs.(ψ̆)),maximum(ψ̆Lims)),max(maximum(abs.(ψ̆)),maximum(ψ̆Lims))]
end 

Colorbar(fig[1,1][1,2],limits=ψ̆Lims,colormap=:bwr,flipaxis=true)


for t=0:100
    @unpack matrices = load(datadir(folderName,"frames","matrices$(@sprintf("%03d", t)).jld2"))
    @unpack params = load(datadir(folderName,"frames","params$(@sprintf("%03d", t)).jld2"))
    @unpack B, Bᵀ, C, R, cellPositions = matrices
    @unpack nCells = params

    cellPolygons = makeCellPolygons(params,matrices)
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
