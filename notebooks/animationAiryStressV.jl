# Script to produce a movie of Airy Stress evolution over vertex network for a given system

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

folderName = "nCells=751_pressureExternal=0.5_realTimetMax=173000.0_stiffnessFactor=10.0_24-03-04-10-11-13"

potentials = Vector{Float64}[]
cellPolygonVectors = Vector{Vector{Point{2,Float64}}}[]
linkTriangleVectors = Vector{Vector{Point{2,Float64}}}[]
files = [datadir("sims", folderName, "frameData", f) for f in readdir(datadir("sims", folderName, "frameData")) if occursin(".jld2", f)]
for t = 2:length(files)
    @show t
    @unpack R, matrices, params = load(files[t];
        typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer" => ParametersContainer,
            "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer" => MatricesContainer))
    ψ̆, spectrum = psivPotential(R, params, matrices)
    push!(potentials, ψ̆)
    cellPolygons = makeCellPolygonsOld(R, params, matrices)
    push!(cellPolygonVectors, cellPolygons)
    linkTriangles = makeLinkTriangles(R, params, matrices)
    push!(linkTriangleVectors, linkTriangles)
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
    @show t
    empty!(ax)
    for k = 1:length(linkTriangleVectors[t])
        poly!(ax,
            linkTriangleVectors[t][k],
            color=potentials[t][k],
            colorrange=ψ̆Lims,
            colormap=:bwr,
            strokewidth=0,
            strokecolor=(:black, 0.0),
        )
    end
    reset_limits!(ax)
    recordframe!(mov)
end

save(datadir("sims", folderName, "movieAiryStressV.mp4"), mov)
