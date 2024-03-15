
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

# MCCtensions = mean(matrices.cellTensions[MCCs])
MCCtensions = matrices.cellTensions[MCCs]
@show MCCtensions
# MCCneighbourTensions = mean(matrices.cellTensions[MCCneighbours])
MCCneighbourTensions = matrices.cellTensions[MCCneighbours]
@show MCCneighbourTensions
# otherTensions = mean(matrices.cellTensions[Not([MCCneighbours...,MCCs...])])
otherTensions = matrices.cellTensions[Not([MCCneighbours...,MCCs...])]
@show otherTensions

minTension = minimum([MCCtensions...,MCCneighbourTensions...,otherTensions...])
maxTension = maximum([MCCtensions...,MCCneighbourTensions...,otherTensions...])

#%%

fig = Figure(size=(500,1000))
# ax1 = Axis(fig[1,1])
# hist1 = fit(Histogram, MCCtensions,  0.0:0.005:2.1)
# barplot!(ax1, hist1, label="MCCs", color=(:red,0.8))
# xlims!(ax1,(1.85,2.1))
# Label(fig[1,1,Bottom()], "MCCs", fontsize = 16)
ax2 = Axis(fig[1,1])
hist2 = fit(Histogram, MCCneighbourTensions,  0.0:0.005:0.25)
barplot!(ax2, hist2, label="Neighbours", color=(:green,0.8))
xlims!(ax2,(0.0,0.25))
Label(fig[1,1,Bottom()], "Neighbours", fontsize = 16)
ax3 = Axis(fig[2,1])
hist3 = fit(Histogram, otherTensions,  0.0:0.005:0.25)
barplot!(ax3, hist3, label="Others", color=(:blue,0.8))
xlims!(ax3,(0.0,0.25))
Label(fig[2,1,Bottom()], "Others", fontsize = 16)


display(fig)
save("Histograms.png",fig)