
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
using StatsBase 

@from "$(projectdir())/src/VertexModelContainers.jl" using VertexModelContainers
@from "$(projectdir())/src/AnalysisFunctions.jl" using AnalysisFunctions

folderName = "nCells=751_pressureExternal=0.5_realTimetMax=173000.0_stiffnessFactor=10.0_24-03-12-15-24-12"

files = [datadir("sims", folderName, "frameData", f) for f in readdir(datadir("sims", folderName, "frameData")) if occursin(".jld2",f)]
@unpack R, matrices, params = load(files[end]; 
    typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer"=>ParametersContainer, 
    "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer"=>MatricesContainer))


cellNeighbourMatrix = matrices.B*matrices.B'
MCCs = findall(x->x>1.5, matrices.Γ)
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

fig = Figure(size=(1000,1000))
ax1 = Axis(fig[1,1])
hist1 = fit(Histogram, nbins=100, matrices.cellTensions[MCCneighbours])
barplot!(ax1, hist1, label="Neighbour cell\ntensions", color=(:green,0.5))
hist2 = fit(Histogram, nbins=100, matrices.cellTensions[Not([MCCneighbours...,MCCs...])])
barplot!(ax1, hist2, label="Other cell\ntensions", color=(:blue,0.5))
axislegend(ax1, merge = true, unique = true)
save(datadir("sims", folderName, "tensionHistograms.png"), fig)

#%%

fig2 = Figure(size=(1000,1000))
ax1 = Axis(fig2[1,1])
hist1 = fit(Histogram, nbins=100, effectivePressure[MCCneighbours])
barplot!(ax1, hist1, label="Neighbour cell\neffective pressures", color=(:green,0.5))
hist2 = fit(Histogram, nbins=100, effectivePressure[Not([MCCneighbours...,MCCs...])])
barplot!(ax1, hist2, label="Other cell\neffective pressures", color=(:blue,0.5))
axislegend(ax1, merge = true, unique = true)
save(datadir("sims", folderName, "peffHistograms.png"), fig2)

#%%

fig3 = Figure(size=(1000,1000))
ax1 = Axis(fig3[1,1])
hist1 = fit(Histogram, nbins=100, shrs[MCCneighbours])
barplot!(ax1, hist1, label="Neighbour cell\nshears", color=(:green,0.5))
hist2 = fit(Histogram, nbins=100, shrs[Not([MCCneighbours...,MCCs...])])
barplot!(ax1, hist2, label="Other cell\nshears", color=(:blue,0.5))
axislegend(ax1, merge = true, unique = true)
save(datadir("sims", folderName, "shearHistograms.png"), fig3)