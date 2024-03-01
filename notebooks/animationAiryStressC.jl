# Script to produce a movie of Airy Stress evolution over cell network for a given system

using JLD2
using SparseArrays
using LinearAlgebra
using DrWatson
using FromFile
using UnPack
using CairoMakie
using Printf
using Colors

@from "$(projectdir())/src/VertexModelContainers.jl" using VertexModelContainers
@from "$(projectdir())/src/AnalysisFunctions.jl" using AnalysisFunctions
@from "$(projectdir())/src/Potentials.jl" using Potentials

# folderName = "/Users/christopher/Postdoc/Code/VertexModel/data/sims/nCells=400_pressureExternal=0.1_realTimetMax=173000.0_stiffnessFactor=5.0_24-02-28-16-00-39"

potentials = Vector{Float64}[]
cellPolygonVectors = Vector{Vector{Point{2,Float64}}}[]
stiffnesses = Vector{Float64}[]
for t = 2:100
    @unpack R, matrices, params = load(datadir(folderName, "frameData", "systemData$(@sprintf("%03d", t)).jld2"))
    ψ̆, spectrum = psicPotential(R, params, matrices)
    push!(potentials, ψ̆)
    cellPolygons = makeCellPolygonsOld(R, params, matrices)
    push!(cellPolygonVectors, cellPolygons)
    push!(stiffnesses, matrices.μ)
end

fig = CairoMakie.Figure(size=(1000, 1000))
ax = Axis(fig[1, 1][1, 1], aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
mov = VideoStream(fig, framerate=5)
globalMax = maximum([maximum(abs.(x)) for x in potentials])
ψ̆Lims = [-globalMax, globalMax]
Colorbar(fig[1, 1][1, 2], limits=ψ̆Lims, colormap=:bwr, flipaxis=true)

for t = 1:length(potentials)
    empty!(ax)
    for i = 1:length(cellPolygonVectors[t])
        if stiffnesses[t][i] < 1.5
            poly!(ax, cellPolygonVectors[t][i], color=potentials[t][i], colormap=:bwr, colorrange=Tuple(ψ̆Lims), strokecolor=(:black, 0.0), strokewidth=0)
        end
    end
    for i = 1:length(cellPolygonVectors[t])
        if stiffnesses[t][i] > 1.5
            poly!(ax, cellPolygonVectors[t][i], color=potentials[t][i], colormap=:bwr, colorrange=Tuple(ψ̆Lims), strokecolor=(:black, 1.0), strokewidth=2)
        end
    end
    reset_limits!(ax)
    recordframe!(mov)
end

save(datadir(folderName, "movieAiryStressC.mp4"), mov)
