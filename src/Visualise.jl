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

        # @show (getRandomColor(i))
        mesh!(ax, connectedVerts, connectedFaces, color=RGB(rand(Xoshiro(i),3)...), shading=NoShading)
        # poly!(ax, connectedVerts, connectedFaces, color=RGB(rand(Xoshiro(i),3)...))
        # poly!(ax, connectedVerts, connectedFaces; color=(:black, 0.5), strokecolor=(:black, 1.0), strokewidth=2)

    end    

    # for i=1:nCells
    #     lines!(Point{3,Float64}.(R[cellVertexOrders[i][1:cellEdgeCount[i]+1]]))
    # end

    # Set limits
    reset_limits!(ax)
    zlims!(ax, (-0.5,0.5))

    recordframe!(mov)

    return nothing

end

export visualise

end
