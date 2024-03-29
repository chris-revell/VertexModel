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

folderName = "nCells=751_pressureExternal=0.5_realTimetMax=173000.0_stiffnessFactor=10.0_24-03-04-10-11-13"

effectivePressureVectors = Vector{Float64}[]
notExcludedCellVectors = Vector{Bool}[]
cellPolygonVectors = Vector{Vector{Point{2,Float64}}}[]
files = [datadir("sims", folderName, "frameData", f) for f in readdir(datadir("sims", folderName, "frameData")) if occursin(".jld2",f)]
for t = 2:length(files)
    @show t
    @unpack R, matrices, params = load(files[t]; 
        typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer"=>ParametersContainer, 
        "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer"=>MatricesContainer))
    effectivePressure = effectiveCellPressure(matrices.cellPressures, matrices.cellTensions, matrices.cellPerimeters, matrices.cellAreas)
    push!(effectivePressureVectors, effectivePressure)
    notExcludedCells = fill(true, params.nCells)
    for j in findall(x -> x != 0, matrices.boundaryEdges)
        notExcludedCells[findnz(matrices.B[:, j])[1][1]] = false
    end
    push!(notExcludedCellVectors, notExcludedCells)
    cellPolygons = makeCellPolygonsOld(R, params, matrices)
    push!(cellPolygonVectors, cellPolygons)
end

fig = CairoMakie.Figure(size=(1000, 1000))
ax = Axis(fig[1, 1][1, 1], aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
mov = VideoStream(fig, framerate=5)
globalPeffMin = minimum([minimum(effectivePressureVectors[t][notExcludedCellVectors[t]]) for t = 1:length(effectivePressureVectors)])
globalPeffMax = maximum([maximum(effectivePressureVectors[t][notExcludedCellVectors[t]]) for t = 1:length(effectivePressureVectors)])
globalLimit = max(abs(globalPeffMin), abs(globalPeffMax))
pLims = (-globalLimit, globalLimit)
Colorbar(fig[1, 1][1, 2], colormap=:bwr, limits=pLims, flipaxis=true)

for t = 1:length(effectivePressureVectors)
    empty!(ax)
    for i = 1:length(effectivePressureVectors[t])
        if notExcludedCellVectors[t][i]
            poly!(ax,
                cellPolygonVectors[t][i],
                color=effectivePressureVectors[t][i],
                colormap=:bwr,
                colorrange=pLims,
                strokecolor=(:black,0.0),
                strokewidth=0,
            )
        end
    end
    reset_limits!(ax)
    recordframe!(mov)
end

save(datadir("sims", folderName, "movieEffectivePressures.mp4"), mov)
