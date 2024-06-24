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

folderName = "pressureExternal=0.5_stiffnessFactor=2.0_γ=0.2_24-06-18-17-39-57"

cellPressureVectors = Vector{Float64}[]
notExcludedCellVectors = Vector{Bool}[]
MCCs = Vector{Int64}[]
cellPolygonVectors = Vector{Vector{Point{2,Float64}}}[]
files = [datadir("sims", folderName, "frameData", f) for f in readdir(datadir("sims", folderName, "frameData")) if occursin(".jld2",f)]
for t = 5:length(files)
    @show t
    @unpack R, matrices, params = load(files[t]; 
        typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer"=>ParametersContainer, 
        "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer"=>MatricesContainer))
    # effectivePressure = effectiveCellPressure(matrices.cellPressures, matrices.cellTensions, matrices.cellPerimeters, matrices.cellAreas)
    push!(cellPressureVectors, matrices.cellPressures)
    notExcludedCells = fill(true, params.nCells)
    for j in findall(x -> x != 0, matrices.boundaryEdges)
        notExcludedCells[findnz(matrices.B[:, j])[1][1]] = false
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
# globalPressuresMin = minimum([minimum(cellPressureVectors[t][notExcludedCellVectors[t]]) for t = 1:length(cellPressureVectors)])
globalPressuresMin = minimum([minimum(cellPressureVectors[t]) for t = 1:length(cellPressureVectors)])
# globalPressuresMax = maximum([maximum(cellPressureVectors[t][notExcludedCellVectors[t]]) for t = 1:length(cellPressureVectors)])
globalPressuresMax = maximum([maximum(cellPressureVectors[t]) for t = 1:length(cellPressureVectors)])
# globalLimit = max(abs(globalPeffMin), abs(globalPeffMax))
# pLims = (-globalLimit, globalLimit)
pLims = (globalPressuresMin, globalPressuresMax)
Colorbar(fig[1, 1][1, 2], colormap=:batlow, limits=pLims, flipaxis=true)

for t = 1:length(cellPressureVectors)
    empty!(ax)
    for i = 1:length(cellPressureVectors[t])
        # if notExcludedCellVectors[t][i]            
        if MCCs[t][i] == 0
            poly!(ax,
                cellPolygonVectors[t][i],
                color=cellPressureVectors[t][i],
                colormap=:batlow,
                colorrange=pLims,
                strokecolor=(:black,0.0),
                strokewidth=2)
        else 
            poly!(ax,
                cellPolygonVectors[t][i],
                color=cellPressureVectors[t][i],
                colormap=:batlow,
                colorrange=pLims,
                strokecolor=(:black,1.0),
                strokewidth=3)
        end
        # else
        #     poly!(ax,
        #         cellPolygonVectors[t][i],
        #         color=(:black,0.5),
        #         strokecolor=(:black,0.0),
        #         strokewidth=2)
        # end
    end
    reset_limits!(ax)
    recordframe!(mov)
    # t==72 ? save(datadir("sims", folderName, "cellPressures072.png"), fig) : nothing 
end

# save(datadir("sims", folderName, "cellPressures100.png"), fig)
save(datadir("sims", folderName, "movieCellPressures.mp4"), mov)
