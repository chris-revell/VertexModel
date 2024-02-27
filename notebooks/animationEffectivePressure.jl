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
folderName = "/Users/christopher/Postdoc/Code/VertexModel/data/sims/L₀=0.75_nCells=800_realTimetMax=17300.0_γ=0.2_24-02-27-08-48-35"

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
    # notExcludedCells[matrices.μ.>1.5] .= false
    for j in findall(x->x!=0, matrices.boundaryEdges)
        notExcludedCells[findnz(matrices.B[:,j])[1][1]] = false
    end
    push!(notExcludedCellVectors,notExcludedCells)    
end 

globalPeffMin = minimum([minimum(effectivePressureVectors[t][notExcludedCellVectors[t]]) for t=1:101])
globalPeffMax = maximum([maximum(effectivePressureVectors[t][notExcludedCellVectors[t]]) for t=1:101])
globalLimit = max(abs(globalPeffMin),abs(globalPeffMax))
pLims = [-globalLimit,globalLimit]

Colorbar(fig[1,1][1,2],colormap=:bwr,limits=pLims,flipaxis=true)

for t=0:100
    @unpack R, matrices, params = load(datadir(folderName,"frameData","systemData$(@sprintf("%03d", t)).jld2"))
    @unpack nCells = params

    cellPolygons = makeCellPolygonsOld(R,params,matrices)
    empty!(ax)
    for i=1:nCells
        if notExcludedCellVectors[t+1][i]
            if matrices.μ[i]>1.5
                # poly!(ax,cellPolygons[i],color=effectivePressureVectors[t+1][i],colormap=:bwr,colorrange=Tuple(pLims), strokecolor=(:black,1.0),strokewidth=2)
            else 
                poly!(ax,cellPolygons[i],color=effectivePressureVectors[t+1][i],colormap=:bwr,colorrange=Tuple(pLims), strokecolor=(:black,0.2),strokewidth=2)
            end
        else
            poly!(ax,cellPolygons[i],color=(:black,0.5), strokecolor=(:black,0.2),strokewidth=2)
        end
    end
    for i=1:nCells
        if matrices.μ[i]>1.5
            poly!(ax,cellPolygons[i],color=effectivePressureVectors[t+1][i],colormap=:bwr,colorrange=Tuple(pLims), strokecolor=(:black,1.0),strokewidth=2)
        end
    end
    reset_limits!(ax)
    recordframe!(mov)
end

save(datadir(folderName,"movieEffectivePressures.mp4"),mov)

# end