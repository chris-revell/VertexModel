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

folderName = "pressureExternal=0.5_stiffnessFactor=10.0_γ=0.2_24-06-24-16-05-11"

effectivePressureVectors = Vector{Float64}[]
notExcludedCellVectors = Vector{Bool}[]
MCCs = Vector{Int64}[]
cellPolygonVectors = Vector{Vector{Point{2,Float64}}}[]
files = [datadir("sims", folderName, "frameData", f) for f in readdir(datadir("sims", folderName, "frameData")) if occursin(".jld2",f)]
for t = 5:length(files)
    # @show t
    @unpack R, matrices, params = load(files[t]; 
        typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer"=>ParametersContainer, 
        "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer"=>MatricesContainer))
    effectivePressure = effectiveCellPressure(matrices.cellPressures, matrices.cellTensions, matrices.cellPerimeters, matrices.cellAreas)
    push!(effectivePressureVectors, effectivePressure)
    notExcludedCells = fill(true, params.nCells)
    for j in findall(x -> x != 0, matrices.boundaryEdges)
        notExcludedCells[findall(x->x!=0, matrices.B[:, j])[1]] = false
    end
    for k in findall(x -> x != 0, matrices.MCCsList)
        notExcludedCells[k] = false
    end
    push!(notExcludedCellVectors, notExcludedCells)
    cellPolygons = makeCellPolygonsOld(R, params, matrices)
    push!(cellPolygonVectors, cellPolygons)
    push!(MCCs, matrices.MCCsList)
end

fig = CairoMakie.Figure(size=(1000, 1000))
ax = Axis(fig[1, 1][1, 1], aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
mov = VideoStream(fig, framerate=5)
globalPeffMin = minimum([minimum(effectivePressureVectors[t][notExcludedCellVectors[t]]) for t = 1:length(effectivePressureVectors)])
# globalPeffMin = minimum([minimum(effectivePressureVectors[t]) for t = 1:length(effectivePressureVectors)])
globalPeffMax = maximum([maximum(effectivePressureVectors[t][notExcludedCellVectors[t]]) for t = 1:length(effectivePressureVectors)])
# globalPeffMax = maximum([maximum(effectivePressureVectors[t]) for t = 1:length(effectivePressureVectors)])

lims = (globalPeffMin, globalPeffMax)

Colorbar(fig[1, 1][1, 2], colormap=:batlow, limits=lims, flipaxis=true)

for t = 1:length(effectivePressureVectors)
    empty!(ax)
    for i = 1:length(effectivePressureVectors[t])
        if notExcludedCellVectors[t][i]
            poly!(ax,
                cellPolygonVectors[t][i],
                color=effectivePressureVectors[t][i],
                colormap=:batlow,
                colorrange=lims,
                strokecolor=(:black,0.5),
                strokewidth=1,
            )
        end
    end
    for i = 1:length(effectivePressureVectors[t])
        if !notExcludedCellVectors[t][i]
            poly!(ax,
                cellPolygonVectors[t][i],
                color=(:black,0.5),
                strokecolor=(:black,(MCCs[t][i] != 0 ? 1.0 : 0.0)),
                strokewidth=3,
            )
        end
    end
    reset_limits!(ax)
    recordframe!(mov)
    # t==72 ? save(datadir("sims", folderName, "cellEffectivePressures072.png"), fig) : nothing 
end
# save(datadir("sims", folderName, "cellEffectivePressures100.png"), fig)
save(datadir("sims", folderName, "movieCellEffectivePressures.mp4"), mov)


