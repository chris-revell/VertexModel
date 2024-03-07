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

folderName = "/Users/christopher/Postdoc/Code/VertexModel/data/sims/nCells=751_pressureExternal=0.5_realTimetMax=173000.0_stiffnessFactor=10.0_24-03-04-10-11-13"# folderName = "/Users/christopher/Postdoc/Code/VertexModel/data/sims/nCells=400_pressureExternal=0.1_realTimetMax=173000.0_stiffnessFactor=5.0_24-02-28-16-00-39"

potentials = Vector{Float64}[]
notExcludedVertVectors = Vector{Bool}[]
cellPolygonVectors = Vector{Vector{Point{2,Float64}}}[]
linkTriangleVectors = Vector{Vector{Point{2,Float64}}}[]
files = [datadir(folderName, "frameData", f) for f in readdir(datadir(folderName, "frameData")) if occursin(".jld2",f)]
for t = 2:length(files)
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
    reset_limits!(ax)
    recordframe!(mov)
end

save(datadir(folderName, "movieMindlinStressV.mp4"), mov)
