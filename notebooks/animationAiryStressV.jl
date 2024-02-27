# Script to produce a movie of Airy Stress evolution over vertex network for a given system

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
@from "$(projectdir())/src/AnalysisFunctions.jl" using AnalysisFunctions
@from "$(projectdir())/src/Laplacians.jl" using Laplacians
@from "$(projectdir())/src/Potentials.jl" using Potentials

# for f in [f for f in readdir(datadir("sims/examples")) if occursin("γ",f)]

folderName = "sims/L₀=0.75_nCells=400_realTimetMax=43200.0_γ=0.2_24-02-21-10-24-23"

fig = CairoMakie.Figure(size=(1000,1000))
ax = Axis(fig[1,1][1,1],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
mov = VideoStream(fig, framerate=5)

potentials = Vector{Float64}[]
for t=0:100
    @unpack R, matrices, params = load(datadir(folderName,"frameData","systemData$(@sprintf("%03d", t)).jld2"))
    ψ̆, spectrum = psivPotential(R,params,matrices)
    push!(potentials,ψ̆)
end 

globalMax = maximum([maximum(abs.(x)) for x in potentials])
ψ̆Lims = [-globalMax,globalMax]

Colorbar(fig[1,1][1,2],limits=ψ̆Lims,colormap=:bwr,flipaxis=true)

for t=0:100
    @unpack R, matrices, params = load(datadir(folderName,"frameData","systemData$(@sprintf("%03d", t)).jld2"))
    @unpack nCells, nVerts = params
    
    cellPolygons = makeCellPolygonsOld(R,params,matrices)
    linkTriangles = makeLinkTriangles(R,params,matrices)
    empty!(ax)
    # for k=1:nVerts
    #     poly!(ax,linkTriangles[k],color=potentials[t+1][k],colorrange=ψ̆Lims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25))
    # end
    for i=1:nCells
        if matrices.μ[i]>1.5
            # poly!(ax,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=2)
        else
            poly!(ax,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,0.2),strokewidth=2)
        end
    end
    for i=1:nCells
        if matrices.μ[i]>1.5
            poly!(ax,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=2)
        end
    end
    reset_limits!(ax)
    recordframe!(mov)
end

save(datadir(folderName,"movieAiryStressV.mp4"),mov)

# end