
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

# folderName = "MCCComparison/pressureExternal=0.5_stiffnessFactor=10.0_γ=0.2_24-06-24-21-29-00"

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

fig = Figure(size=(1000,1000))
ax = Axis(fig[1,1])
# tensionbins = minimum(matrices.cellTensions):0.1:maximum(matrices.cellTensions)
tensionbins = 0.0:0.01:1.0 #maximum(matrices.cellTensions)
hist1 = fit(Histogram, matrices.cellTensions[MCCneighbours], tensionbins, closed=:right)
barplot!(ax, hist1, label="Neighbour cell\ntensions", color=(:green,0.5))
hist2 = fit(Histogram, matrices.cellTensions[Not([MCCneighbours...,MCCs...])], tensionbins, closed=:right)
barplot!(ax, hist2, label="Other cell\ntensions", color=(:blue,0.5))
# axislegend(ax)#, merge = true, unique = true)
save(datadir("sims", folderName, "tensionHistograms.png"), fig)
fig = Figure(size=(1000,1000)); ax = Axis(fig[1,1])
hist1Normalised = normalize(hist1, mode=:pdf)
hist2Normalised = normalize(hist2, mode=:pdf)
barplot!(ax, hist1Normalised, label="Neighbour cell\ntensions", color=(:green,0.5))
barplot!(ax, hist2Normalised, label="Other cell\ntensions", color=(:blue,0.5))
# axislegend(ax)#, merge = true, unique = true)
save(datadir("sims", folderName, "tensionDensities.png"), fig)

fig = Figure(size=(1000,1000)); ax = Axis(fig[1,1])
for i=2:length(hist1Normalised.weights)
    hist1Normalised.weights[i]=hist1Normalised.weights[i]+hist1Normalised.weights[i-1]
end
for i=2:length(hist2Normalised.weights)
    hist2Normalised.weights[i]=hist2Normalised.weights[i]+hist2Normalised.weights[i-1]
end
barplot!(ax, hist1Normalised, label="Neighbour cell\ntensions", color=(:green,0.5))
barplot!(ax, hist2Normalised, label="Other cell\ntensions", color=(:blue,0.5))
# axislegend(ax)#, merge = true, unique = true)
save(datadir("sims", folderName, "tensionCumulativeDensities.png"), fig)

#%%

fig = Figure(size=(1000,1000)); ax = Axis(fig[1,1])
hist1 = fit(Histogram, nbins=100, effectivePressure[MCCneighbours], closed=:right)
barplot!(ax, hist1, label="Neighbour cell\neffective pressures", color=(:green,0.5))
hist2 = fit(Histogram, nbins=100, effectivePressure[Not([MCCneighbours...,MCCs...])], closed=:right)
barplot!(ax, hist2, label="Other cell\neffective pressures", color=(:blue,0.5))
axislegend(ax)#, merge = true, unique = true)
save(datadir("sims", folderName, "peffHistograms.png"), fig)
fig = Figure(size=(1000,1000)); ax = Axis(fig[1,1])
hist1Normalised = normalize(hist1, mode=:pdf)
hist2Normalised = normalize(hist2, mode=:pdf)
barplot!(ax, hist1Normalised, label="Neighbour cell\neffective pressures", color=(:green,0.5))
barplot!(ax, hist2Normalised, label="Other cell\neffective pressures", color=(:blue,0.5))
axislegend(ax)#, merge = true, unique = true)
save(datadir("sims", folderName, "peffDensities.png"), fig)

fig = Figure(size=(1000,1000)); ax = Axis(fig[1,1])
for i=2:length(hist1Normalised.weights)
    hist1Normalised.weights[i]=hist1Normalised.weights[i]+hist1Normalised.weights[i-1]
end
for i=2:length(hist2Normalised.weights)
    hist2Normalised.weights[i]=hist2Normalised.weights[i]+hist2Normalised.weights[i-1]
end
barplot!(ax, hist1Normalised, label="Neighbour cell\neffective pressures", color=(:green,0.5))
barplot!(ax, hist2Normalised, label="Other cell\neffective pressures", color=(:blue,0.5))
axislegend(ax)#, merge = true, unique = true)
save(datadir("sims", folderName, "peffCumulativeDensities.png"), fig)
#%%
#%%

fig = Figure(size=(1000,1000)); ax = Axis(fig[1,1])
hist1 = fit(Histogram, nbins=100, matrices.cellPressures[MCCneighbours], closed=:right)
barplot!(ax, hist1, label="Neighbour cell\neffective pressures", color=(:green,0.5))
hist2 = fit(Histogram, nbins=100, matrices.cellPressures[Not([MCCneighbours...,MCCs...])], closed=:right)
barplot!(ax, hist2, label="Other cell\neffective pressures", color=(:blue,0.5))
axislegend(ax)#, merge = true, unique = true)
save(datadir("sims", folderName, "peffHistograms.png"), fig)
fig = Figure(size=(1000,1000)); ax = Axis(fig[1,1])
hist1Normalised = normalize(hist1, mode=:pdf)
hist2Normalised = normalize(hist2, mode=:pdf)
barplot!(ax, hist1Normalised, label="Neighbour cell\n pressures", color=(:green,0.5))
barplot!(ax, hist2Normalised, label="Other cell\n pressures", color=(:blue,0.5))
axislegend(ax)#, merge = true, unique = true)
save(datadir("sims", folderName, "pressureDensities.png"), fig)

fig = Figure(size=(1000,1000)); ax = Axis(fig[1,1])
for i=2:length(hist1Normalised.weights)
    hist1Normalised.weights[i]=hist1Normalised.weights[i]+hist1Normalised.weights[i-1]
end
for i=2:length(hist2Normalised.weights)
    hist2Normalised.weights[i]=hist2Normalised.weights[i]+hist2Normalised.weights[i-1]
end
barplot!(ax, hist1Normalised, label="Neighbour cell\n pressures", color=(:green,0.5))
barplot!(ax, hist2Normalised, label="Other cell\n pressures", color=(:blue,0.5))
axislegend(ax)#, merge = true, unique = true)
save(datadir("sims", folderName, "pressureCumulativeDensities.png"), fig)

#%%

fig = Figure(size=(1000,1000)); ax = Axis(fig[1,1])
hist1 = fit(Histogram, nbins=100, shrs[MCCneighbours], closed=:right)
barplot!(ax, hist1, label="Neighbour cell\nshears", color=(:green,0.5))
hist2 = fit(Histogram, nbins=100, shrs[Not([MCCneighbours...,MCCs...])], closed=:right)
barplot!(ax, hist2, label="Other cell\nshears", color=(:blue,0.5))
axislegend(ax)#, merge = true, unique = true)
save(datadir("sims", folderName, "shearHistograms.png"), fig)
fig = Figure(size=(1000,1000)); ax = Axis(fig[1,1])
hist1Normalised = normalize(hist1, mode=:pdf)
hist2Normalised = normalize(hist2, mode=:pdf)
barplot!(ax, hist1Normalised, label="Neighbour cell\nshears", color=(:green,0.5))
barplot!(ax, hist2Normalised, label="Other cell\nshears", color=(:blue,0.5))
axislegend(ax)#, merge = true, unique = true)
save(datadir("sims", folderName, "shearDensities.png"), fig)


fig = Figure(size=(1000,1000)); ax = Axis(fig[1,1])
for i=2:length(hist1Normalised.weights)
    hist1Normalised.weights[i]=hist1Normalised.weights[i]+hist1Normalised.weights[i-1]
end
for i=2:length(hist2Normalised.weights)
    hist2Normalised.weights[i]=hist2Normalised.weights[i]+hist2Normalised.weights[i-1]
end
barplot!(ax, hist1Normalised, label="Neighbour cell\nshears", color=(:green,0.5))
barplot!(ax, hist2Normalised, label="Other cell\nshears", color=(:blue,0.5))
axislegend(ax)#, merge = true, unique = true)
save(datadir("sims", folderName, "shearCumulativeDensities.png"), fig)
#%%

# fig7 = Figure(size=(1000,1000))
# ax7 = Axis(fig7[1,1])
# range = minimum(shrs):0.01:maximum(shrs)
# cDif1 = ecdf(shrs[MCCneighbours])
# lines!(ax7, collect(range), cDif1(range), label="Neighbour cell\nshears", color=(:green,0.5))
# cDif2 = ecdf(shrs[Not([MCCneighbours...,MCCs...])])
# lines!(ax7, collect(range), cDif2(range), label="Other cell\nshears", color=(:blue,0.5))
# axislegend(ax7, merge = true, unique = true)
# display(fig7)
# save(datadir("sims", folderName, "shearECDF.png"), fig7)


# fig8 = Figure(size=(1000,1000))
# ax8 = Axis(fig8[1,1])
# hist1Normalised = normalize(hist1, mode=:pdf)
# hist2Normalised = normalize(hist2, mode=:pdf)
# barplot!(ax6, hist1Normalised, label="Neighbour cell\nshears", color=(:green,0.5))
# barplot!(ax6, hist2Normalised, label="Other cell\nshears", color=(:blue,0.5))
# axislegend(ax6, merge = true, unique = true)
# save(datadir("sims", folderName, "shearDensities.png"), fig6)