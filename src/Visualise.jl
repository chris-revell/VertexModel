#
#  Visualise.jl
#  VertexModel
#
#  Created by Christopher Revell on 16/02/2021.
#
#
#

module Visualise

# Julia packages
using Printf
using LinearAlgebra
using ColorSchemes
using Colors
using UnPack
using GeometryBasics
using Random
using Makie
using GLMakie
using CairoMakie
using StaticArrays
using SparseArrays
using CircularArrays
using FromFile
using DrWatson
using Random
using Distributions

# Local modules
@from "OrderAroundCell.jl" using OrderAroundCell
@from "AnalysisFunctions.jl" using AnalysisFunctions

function visualise(R, t, fig, ax, mov, params, matrices)

    @unpack cellEdgeCount,
        cellVertexOrders,
        cellEdgeOrders,
        cellPositions,
        edgeMidpoints,
        F,
        edgeMidpointLinks,
        μ = matrices
    @unpack nEdges,
        nVerts,
        nCells = params

    empty!(ax)

    # ax.title = "t = $(@sprintf("%.2f", t))"

    for i=1:nCells
        verts = Float64[]
        for k=1:length(cellVertexOrders[i])
            append!(verts, R[cellVertexOrders[i][k]])
            append!(verts, R[cellVertexOrders[i][k+1]])
            append!(verts, cellPositions[i])
        end
        connectedVerts = connect(verts, Point{3})
        connectedFaces = connect(1:length(connectedVerts), TriangleFace)
        mesh!(ax, connectedVerts, connectedFaces, color=RGB(rand(Xoshiro(i),3)...), shading=NoShading)
    end    

    
    # text!(ax, Point{3,Float64}.(matrices.cellPositions), text = string.(collect(1:nCells)), color=:red)
    # scatter!(ax, Point{3,Float64}.(matrices.cellPositions), color=:red)
    # text!(ax, Point{3,Float64}.(matrices.edgeMidpoints), text = string.(collect(1:nEdges)), color=:green)
    # scatter!(ax, Point{3,Float64}.(matrices.edgeMidpoints), color=:green)


    # Set limits
    reset_limits!(ax)

    recordframe!(mov)

    return nothing

end


function visualise3DInstance(R, params, matrices; labels=false)
    @unpack cellEdgeCount,
        cellVertexOrders,
        cellEdgeOrders,
        cellPositions,
        edgeMidpoints,
        F,
        edgeMidpointLinks,
        μ = matrices
    @unpack nEdges,
        nVerts,
        nCells = params
    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    GLMakie.activate!()
    fig = GLMakie.Figure(size=(1000,1000))
    ax = LScene(fig[1,1])

    for i=1:nCells
        verts = Float64[]
        for k=1:length(cellVertexOrders[i])
            append!(verts, R[cellVertexOrders[i][k]])
            append!(verts, R[cellVertexOrders[i][k+1]])
            append!(verts, cellPositions[i])
        end
        connectedVerts = connect(verts, Point{3})
        connectedFaces = connect(1:length(connectedVerts), TriangleFace)
        mesh!(ax, connectedVerts, connectedFaces, color=RGB(rand(Xoshiro(i),3)...), shading=NoShading)
    end    
    if labels 
        text!(ax, Point{3,Float64}.([p.+[0.0,0.0,0.1] for p in matrices.cellPositions]), text = string.(collect(1:nCells)), color=:red)
        scatter!(ax, Point{3,Float64}.([p.+[0.0,0.0,0.1] for p in matrices.cellPositions]), color=:red)
        text!(ax, Point{3,Float64}.([p.+[0.0,0.0,0.1] for p in matrices.edgeMidpoints]), text = string.(collect(1:nEdges)), color=:green)
        scatter!(ax, Point{3,Float64}.([p.+[0.0,0.0,0.1] for p in matrices.edgeMidpoints]), color=:green)
        text!(ax, Point{3,Float64}.([p.+[0.0,0.0,0.1] for p in R]), text = string.(collect(1:nVerts)), color=:blue)
        scatter!(ax, Point{3,Float64}.([p.+[0.0,0.0,0.1] for p in R]), color=:blue)
    end
    reset_limits!(ax)
    return fig
end


function visualiseSet(subFolder)

    files = [f for f in readdir(datadir("sims", subFolder, "frameData")) if occursin(".jld2", f)]
    all_R = []
    all_matrices = []
    all_params = []
    for f in files
        data = load(datadir("sims", subFolder, "frameData", f))
        push!(all_R, data["R"])
        push!(all_matrices, data["matrices"])
        push!(all_params, data["params"])
    end

    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    fig = Figure(size=(1000,1000))
    ax = Axis3(fig[1,1], aspect=:data)
    mov = VideoStream(fig, framerate=5)
    for t=1:length(all_R)
        empty!(ax)
        for i=1:all_params[t].nCells
            verts = Float64[]
            for k=1:length(all_matrices[t].cellVertexOrders[i])
                append!(verts, all_R[t][all_matrices[t].cellVertexOrders[i][k]])
                append!(verts, all_R[t][all_matrices[t].cellVertexOrders[i][k+1]])
                append!(verts, all_matrices[t].cellPositions[i])
            end
            connectedVerts = connect(verts, Point{3})
            connectedFaces = connect(1:length(connectedVerts), TriangleFace)
            mesh!(ax, connectedVerts, connectedFaces, color=RGB(rand(Xoshiro(i),3)...), shading=NoShading)
        end   
        reset_limits!(ax)
        recordframe!(mov)
    end
    save(datadir("sims", subFolder, "$(splitpath(subFolder)[end]).mp4"), mov)
    return nothing
end

export visualise
export visualise3DInstance
export visualiseSet 

end
