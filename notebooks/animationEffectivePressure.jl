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

# for f in [f for f in readdir(datadir("sims/examples")) if occursin("γ",f)]

# folderName = "sims/examples/$f"
folderName = "sims/L₀=0.75_nCells=400_realTimetMax=43200.0_γ=0.2_24-02-21-10-24-23"

fig = CairoMakie.Figure(size=(1000,1000))
ax = Axis(fig[1,1][1,1],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
mov = VideoStream(fig, framerate=5)

effectivePressureVectors = Vector{Float64}[]
notExcludedCellVectors = Vector{Bool}[]
for t=0:100
    @unpack R, matrices, params = load(datadir(folderName,"frameData","systemData$(@sprintf("%03d", t)).jld2"))
    effectivePressure = effectiveCellPressure(matrices.cellPressures,matrices.cellTensions,matrices.cellPerimeters,matrices.cellAreas)
    push!(effectivePressureVectors,effectivePressure)
    notExcludedCells = fill(true,params.nCells)
    notExcludedCells[matrices.Γ.>1.5] .= false
    for j in findall(x->x!=0, matrices.boundaryEdges)
        notExcludedCells[findnz(matrices.B[:,j])[1][1]] = false
    end
    push!(notExcludedCellVectors,notExcludedCells)    
end 

globalPeffMin = minimum([minimum(effectivePressureVectors[t][notExcludedCellVectors[t]]) for t=1:101])
globalPeffMax = maximum([maximum(effectivePressureVectors[t][notExcludedCellVectors[t]]) for t=1:101])
pLims = [globalPeffMin,globalPeffMax]

Colorbar(fig[1,1][1,2],colormap=:imola,limits=pLims,flipaxis=true)

for t=0:100
    @unpack R, matrices, params = load(datadir(folderName,"frameData","systemData$(@sprintf("%03d", t)).jld2"))
    @unpack nCells = params

    cellPolygons = makeCellPolygonsOld(R,params,matrices)
    empty!(ax)
    for i=1:nCells
        if notExcludedCellVectors[t+1][i]
            poly!(ax,cellPolygons[i],color=effectivePressureVectorss[t+1][i],colormap=:imola,colorrange=Tuple(pLims), strokecolor=(:black,1.0),strokewidth=2)
        else
            poly!(ax,cellPolygons[i],color=:black, strokecolor=(:black,1.0),strokewidth=2)
        end
    end
    reset_limits!(ax)
    recordframe!(mov)
end

save(datadir(folderName,"movieEffectivePressures.mp4"),mov)

# end