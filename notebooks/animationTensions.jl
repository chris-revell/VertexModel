# Script to produce a movie of Airy Stress evolution over cell network for a given system

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

@from "$(projectdir())/src/VertexModelContainers.jl" using VertexModelContainers
@from "$(projectdir())/src/AnalysisFunctions.jl" using AnalysisFunctions

folderName = "nCells=751_pressureExternal=0.5_realTimetMax=173000.0_stiffnessFactor=10.0_24-03-04-10-11-13"

tensionVectors = Vector{Float64}[]
edgeMidpointVectors = Vector{StaticVector{2,Float64}}[]
edgeTangentVectors = Vector{StaticVector{2,Float64}}[]
notExcludedEdgeVectors = Vector{Bool}[]
cellPolygonVectors = Vector{Vector{Point{2,Float64}}}[]
files = [datadir("sims", folderName, "frameData", f) for f in readdir(datadir("sims", folderName, "frameData")) if occursin(".jld2",f)]
for t = 2:length(files)
    @show t
    @unpack R, matrices, params = load(files[t]; 
        typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer"=>ParametersContainer, 
        "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer"=>MatricesContainer))
    
    edgeTensions = matrices.B̄ᵀ*matrices.cellTensions
    
    push!(tensionVectors, edgeTensions)
    
    notExcludedEdges = fill(true, params.nEdges)
    for i in findall(x -> x > 1.5, matrices.μ)
        notExcludedEdges[matrices.cellEdgeOrders[i]] .= false
    end
    push!(notExcludedEdgeVectors, notExcludedEdges)
    cellPolygons = makeCellPolygonsOld(R, params, matrices)
    push!(cellPolygonVectors, cellPolygons)
    push!(edgeMidpointVectors, matrices.edgeMidpoints)
    push!(edgeTangentVectors, matrices.edgeTangents)
end

fig = CairoMakie.Figure(size=(1000, 1000))
ax = Axis(fig[1, 1][1, 1], aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
mov = VideoStream(fig, framerate=5)
globalTensionMin = minimum([minimum(tensionVectors[t][notExcludedEdgeVectors[t]]) for t = 1:length(tensionVectors)])
globalTensionMax = maximum([maximum(tensionVectors[t][notExcludedEdgeVectors[t]]) for t = 1:length(tensionVectors)])
# globalLimit = max(abs(globalTensionMin), abs(globalTensionMax))
# pLims = (-globalLimit, globalLimit)
pLims = (globalTensionMin, globalTensionMax)
Colorbar(fig[1, 1][1, 2], colormap=:hawaii, limits=pLims, flipaxis=true)

for t = 1:length(tensionVectors)
    empty!(ax)
    p1s = Point{2,Float64}.(edgeMidpointVectors[t].-edgeTangentVectors[t]./2.0)
    p2s = Point{2,Float64}.(edgeMidpointVectors[t].+edgeTangentVectors[t]./2.0)
    for j = 1:length(tensionVectors[t])
        if notExcludedEdgeVectors[t][j]
            lines!(ax,
                [p1s[j],p2s[j]],
                color=tensionVectors[t][j],
                colormap=:hawaii,
                colorrange=pLims,
                strokewidth=2,
            )
        else 
            lines!(ax,
                [p1s[j],p2s[j]],
                color=(:black, 1.0),
                strokewidth=2,
            )
        end
    end
    reset_limits!(ax)
    recordframe!(mov)
end

save(datadir("sims", folderName, "movieTensions.mp4"), mov)
