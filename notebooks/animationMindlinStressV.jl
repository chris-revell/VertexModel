# Script to produce a movie of Mindlin Stress evolution over vertex network for a given system

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

folderName = "L₀=0.75_nCells=61_pressureExternal=0.5_realTimetMax=86400.0_stiffnessFactor=1.0_γ=0.2_24-05-08-16-38-39"

potentials = Vector{Float64}[]
stiffnesses = Vector{Float64}[]
notExcludedVertVectors = Vector{Bool}[]
cellPolygonVectors = Vector{Vector{Point{2,Float64}}}[]
linkTriangleVectors = Vector{Vector{Point{2,Float64}}}[]
files = [datadir("sims", folderName, "frameData", f) for f in readdir(datadir("sims", folderName, "frameData")) if occursin(".jld2",f)]
for t = 5:length(files)
    @show t
    @unpack R, matrices, params = load(files[t]; 
        typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer"=>ParametersContainer, 
        "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer"=>MatricesContainer))
    ψ̆, spectrum = capitalPsivPotential(R, params, matrices)
    push!(potentials, ψ̆)
    notExcludedVerts = fill(true, params.nVerts)
    for j in findall(x -> x != 0, matrices.boundaryVertices)
        notExcludedVerts[j] = false
    end
    push!(notExcludedVertVectors, notExcludedVerts)
    cellPolygons = makeCellPolygonsOld(R, params, matrices)
    push!(cellPolygonVectors, cellPolygons)
    linkTriangles = makeLinkTriangles(R, params, matrices)
    push!(linkTriangleVectors, linkTriangles)
    push!(stiffnesses, matrices.μ)
end

fig = CairoMakie.Figure(size=(1000, 1000))
ax = Axis(fig[1, 1][1, 1], aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
mov = VideoStream(fig, framerate=5)
globalMax = maximum([maximum(abs.(potentials[t][notExcludedVertVectors[t]])) for t = 1:length(potentials)])
ψ̆Lims = (-globalMax, globalMax)
Colorbar(fig[1, 1][1, 2], limits=ψ̆Lims, colormap=:bwr, flipaxis=true)

for t = 1:length(potentials)
    @show t
    empty!(ax)
    for k = 1:length(linkTriangleVectors[t])
        if notExcludedVertVectors[t][k]
            poly!(ax,
            linkTriangleVectors[t][k],
            color=potentials[t][k],
            colorrange=ψ̆Lims,
            colormap=:bwr,
            strokewidth=0,
            strokecolor=(:black, 0.0))
        else
            poly!(ax,
            linkTriangleVectors[t][k],
            color=(:black,0.5),
            strokewidth=0,
            strokecolor=(:black,0.0))
        end
    end
    for i = 1:length(cellPolygonVectors[t])
        if stiffnesses[t][i] > 1.5
            poly!(ax,
            cellPolygonVectors[t][i],
            color=(:white,0.0),
            strokecolor=(:black, 1.0),
            strokewidth=2)
        end
    end
    reset_limits!(ax)
    recordframe!(mov)
end

save(datadir("sims", folderName, "movieMindlinStressV.mp4"), mov)
