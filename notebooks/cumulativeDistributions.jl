
using JLD2
using SparseArrays
using StaticArrays
using LinearAlgebra
using DrWatson
using FromFile
using UnPack
using CairoMakie
using Printf
using Colors
using Statistics
using InvertedIndices
# using StatsBase 

@from "$(projectdir())/src/VertexModelContainers.jl" using VertexModelContainers
@from "$(projectdir())/src/AnalysisFunctions.jl" using AnalysisFunctions

folderName = "L₀=0.75_nCells=61_pressureExternal=0.5_realTimetMax=432000.0_stiffnessFactor=2.0_γ=0.2_24-06-04-17-04-41"

files = [datadir("sims", folderName, "frameData", f) for f in readdir(datadir("sims", folderName, "frameData")) if occursin(".jld2",f)]
@unpack R, matrices, params = load(files[end]; 
    typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer"=>ParametersContainer, 
    "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer"=>MatricesContainer))

cellNeighbourMatrix = matrices.B*matrices.B'
MCCs = findall(x->x>1.5, matrices.μ)
MCCneighbours = Int64[]
for i in MCCs
    append!(MCCneighbours, [x for x in findall(x->x!=0, cellNeighbourMatrix[i,:]) if x∉MCCs] )
end
unique!(MCCneighbours)
MCCtensions = matrices.cellTensions[MCCs]
MCCneighbourTensions = matrices.cellTensions[MCCneighbours]
otherTensions = matrices.cellTensions[Not([MCCneighbours...,MCCs...])]
Qs = cellQs(matrices.cellPerimeters, matrices.edgeTangents, matrices.B̄)
shrs = cellShears(matrices.cellTensions, matrices.cellPerimeters, matrices.cellAreas, Qs)
effectivePressure = effectiveCellPressure(matrices.cellPressures, matrices.cellTensions, matrices.cellPerimeters, matrices.cellAreas)

#%%

fig = Figure(size=(1000,1000)); ax = Axis(fig[1,1])
sc1 = ecdfplot!(ax, matrices.cellTensions[MCCneighbours]; color=:red, label="Neighbour cell tensions")
sc2 = ecdfplot!(ax, matrices.cellTensions[Not([MCCneighbours...,MCCs...])]; color=:green, label="Other cell tensions")
# axislegend(ax)#, merge = true, unique = true)
display(fig)

save(datadir("sims", folderName, "tensionCumulativeDensities.png"), fig)

#%%

fig = Figure(size=(1000,1000)); ax = Axis(fig[1,1])
ecdfplot!(ax, effectivePressure[MCCneighbours], label="Neighbour cell\neffective pressures", color=:green)
ecdfplot!(ax, effectivePressure[Not([MCCneighbours...,MCCs...])], label="Other cell\neffective pressures", color=:blue)
# axislegend(ax)#, merge = true, unique = true)
save(datadir("sims", folderName, "peffCumulativeDensities.png"), fig)
#%%
#%%

fig = Figure(size=(1000,1000)); ax = Axis(fig[1,1])
ecdfplot!(ax, matrices.cellPressures[MCCneighbours], label="Neighbour cell\n pressures", color=:green)
ecdfplot!(ax, matrices.cellPressures[Not([MCCneighbours...,MCCs...])], label="Other cell\n pressures", color=:blue)
# axislegend(ax)#, merge = true, unique = true)
save(datadir("sims", folderName, "pressureCumulativeDensities.png"), fig)

#%%

fig = Figure(size=(1000,1000)); ax = Axis(fig[1,1])
ecdfplot!(ax, shrs[MCCneighbours], label="Neighbour cell\nshears", color=:green)
ecdfplot!(ax, shrs[Not([MCCneighbours...,MCCs...])], label="Other cell\nshears", color=:blue)
# axislegend(ax)#, merge = true, unique = true)
save(datadir("sims", folderName, "shearCumulativeDensities.png"), fig)
#%%
