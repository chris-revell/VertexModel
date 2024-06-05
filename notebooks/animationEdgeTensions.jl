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

folderName = "L₀=0.75_nCells=61_pressureExternal=0.5_realTimetMax=432000.0_stiffnessFactor=2.0_γ=0.2_24-06-04-17-04-41"

tensionVectors = Vector{Float64}[]
notExcludedEdgeVectors = Vector{Bool}[]
cellPolygonVectors = Vector{Vector{Point{2,Float64}}}[]
edgeTrapeziumVectors = Vector{Vector{Point{2,Float64}}}[]
files = [datadir("sims", folderName, "frameData", f) for f in readdir(datadir("sims", folderName, "frameData")) if occursin(".jld2",f)]
for t = 5:length(files)
    @show t
    @unpack R, matrices, params = load(files[t]; 
        typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer"=>ParametersContainer, 
        "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer"=>MatricesContainer))
    edgeTensions = matrices.B̄ᵀ*matrices.cellTensions
    push!(tensionVectors, edgeTensions)
    notExcludedEdges = fill(true, params.nEdges)
    # for i in findall(x -> x > 1.5, matrices.μ)
    #     notExcludedEdges[matrices.cellEdgeOrders[i]] .= false
    # end
    for i in findall(x->x==1, matrices.boundaryEdges)
        notExcludedEdges[i] = false
    end
    push!(notExcludedEdgeVectors, notExcludedEdges)
    cellPolygons = makeCellPolygons(R, params, matrices)
    push!(cellPolygonVectors, cellPolygons)
    edgeTrapezia = makeEdgeTrapezia(R,params,matrices)
    push!(edgeTrapeziumVectors, edgeTrapezia)
end

fig = CairoMakie.Figure(size=(1000, 1000))
ax = Axis(fig[1, 1][1, 1], aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
mov = VideoStream(fig, framerate=5)
globalCellTensionMin = minimum([minimum(tensionVectors[t][notExcludedEdgeVectors[t]]) for t = 1:length(tensionVectors)])
globalCellTensionMax = maximum([maximum(tensionVectors[t][notExcludedEdgeVectors[t]]) for t = 1:length(tensionVectors)])
# globalLimit = max(abs(globalTensionMin), abs(globalTensionMax))
# pLims = (-globalLimit, globalLimit)
pLims = (globalCellTensionMin, globalCellTensionMax)
Colorbar(fig[1, 1][1, 2], colormap=:batlow, limits=pLims, flipaxis=true)

for t = 1:length(tensionVectors)
    empty!(ax)
    for j = 1:length(tensionVectors[t])
        if notExcludedEdgeVectors[t][j]
            poly!(ax,
                edgeTrapeziumVectors[t][j],
                color=tensionVectors[t][j],
                colorrange=pLims,
                colormap=:batlow,
                strokewidth=1,
                strokecolor=(:black, 0.25)
            )
        else 
            poly!(ax,
                edgeTrapeziumVectors[t][j],
                color=(:black,1.0),
                colorrange=pLims,
                colormap=:batlow,
                strokewidth=1,
                strokecolor=(:black, 0.25)
            )
        end
    end
    for i = 1:length(cellPolygonVectors[t])
        poly!(ax,
            cellPolygonVectors[t][i],
            color=(:white, 0.0),
            strokecolor=(:black, 1.0),
            strokewidth=2,
        )
    end
    reset_limits!(ax)
    recordframe!(mov)
    t==72 ? save(datadir("sims", folderName, "EdgeTensions072.png"), fig) : nothing 
end

save(datadir("sims", folderName, "EdgeTensions100.png"), fig)
save(datadir("sims", folderName, "movieEdgeTensions.mp4"), mov)
