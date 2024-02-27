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
@from "$(projectdir())/src/AnalysisFunctions.jl" using AnalysisFunctions
@from "$(projectdir())/src/Laplacians.jl" using Laplacians
@from "$(projectdir())/src/Potentials.jl" using Potentials
@from "$(projectdir())/src/SpatialData.jl" using SpatialData

# for f in [f for f in readdir(datadir("sims/examples")) if occursin("γ",f)]

folderName = "sims/L₀=0.75_nCells=400_realTimetMax=43200.0_γ=0.2_24-02-21-10-24-23"

fig = CairoMakie.Figure(size=(1000,1000))
ax = Axis(fig[1,1][1,1],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
mov = VideoStream(fig, framerate=5)

potentials = Vector{Float64}[]
for t=0:100
    @unpack R, matrices, params = load(datadir(folderName,"frameData","systemData$(@sprintf("%03d", t)).jld2"))
    ψ̆, spectrum = psicPotential(R,params,matrices)
    push!(potentials,ψ̆)
end 

globalMax = maximum([maximum(abs.(x)) for x in potentials])
ψ̆Lims = [-globalMax,globalMax]

Colorbar(fig[1,1][1,2],limits=ψ̆Lims,colormap=:bwr,flipaxis=true)

for t=0:100
    @unpack R, matrices, params = load(datadir(folderName,"frameData","systemData$(@sprintf("%03d", t)).jld2"))
    @unpack nCells = params

    cellPolygons = makeCellPolygonsOld(R,params,matrices)
    empty!(ax)
    for i=1:nCells
        if matrices.μ[i]>1.5
            # poly!(ax,cellPolygons[i],color=potentials[t+1][i],colormap=:bwr,colorrange=Tuple(ψ̆Lims), strokecolor=(:black,1.0),strokewidth=2)
        else
            poly!(ax,cellPolygons[i],color=potentials[t+1][i],colormap=:bwr,colorrange=Tuple(ψ̆Lims), strokecolor=(:black,0.2),strokewidth=2)
        end
    end
    for i=1:nCells 
        if matrices.μ[i]>1.5
            poly!(ax,cellPolygons[i],color=potentials[t+1][i],colormap=:bwr,colorrange=Tuple(ψ̆Lims), strokecolor=(:black,1.0),strokewidth=2)
        end
    end
    reset_limits!(ax)
    recordframe!(mov)
    # save(datadir(folderName,"airyStress$(t).png"),fig)
end

save(datadir(folderName,"movieAiryStressC.mp4"),mov)
