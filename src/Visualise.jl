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
using CairoMakie
using StaticArrays
using SparseArrays
using CircularArrays
using FromFile
using DrWatson

# Local modules
@from "OrderAroundCell.jl" using OrderAroundCell
@from "AnalysisFunctions.jl" using AnalysisFunctions

function visualise(R, t, fig, ax, mov, params, matrices, plotCells, scatterEdges, scatterVertices, scatterCells, plotForces, plotEdgeMidpointLinks)

    @unpack cellEdgeCount,
        cellVertexOrders,
        cellEdgeOrders,
        cellPositions,
        edgeMidpoints,
        F,
        edgeMidpointLinks,
        μ,
        boundaryCells = matrices
    @unpack initialSystem,
    nEdges, 
    nVerts, 
    nCells, 
    cellsTypeA, 
    cellsTypeB = params


    empty!(ax)

    ax.title = "t = $(@sprintf("%.3f", t))"

    # Plot cells
    if plotCells == 1
        cellPolygons = makeCellPolygons(R, params, matrices)
        for i = 1:nCells

            if initialSystem == "periodic"
                # Check whether it is on the periodic boundary: 
                if boundaryCells[i]==1

                    N_x = 10
                    N_y = 10

                    num_vertices = length(cellPolygons[i])
                    newCellPolygon = zeros(num_vertices, 2)

                    # Initialise a polygon on the other side of the domain: 
                    oppositePolygon1 = zeros(num_vertices, 2)
                    oppositePolygon2 = zeros(num_vertices, 2)
                    oppositePolygon3 = zeros(num_vertices, 2)
                    oppositePolygon4 = zeros(num_vertices, 2)
                    oppositePolygon5 = zeros(num_vertices, 2)
                    oppositePolygon6 = zeros(num_vertices, 2)
                    oppositePolygon7 = zeros(num_vertices, 2)
                    oppositePolygon8 = zeros(num_vertices, 2)

                    # Flags to see which boundaries are being crossed 
                    flag1=0
                    flag2=0

                    for k = 1:num_vertices

                        if norm(cellPolygons[i][k][1]-cellPositions[i][1]) > N_x/2
                            
                            if cellPolygons[i][k][1] > cellPositions[i][1]
                                newCellPolygon[k,1] = cellPolygons[i][k][1] - N_x
                                flag1 = 1
                            else
                                newCellPolygon[k,1] = cellPolygons[i][k][1] + N_x
                                flag1 = 2
                            end
                        else
                            newCellPolygon[k,1] = cellPolygons[i][k][1]
                        end

                        if norm(cellPolygons[i][k][2]-cellPositions[i][2]) > N_y/2
                            
                            if cellPolygons[i][k][2] > cellPositions[i][2]
                                flag2 = 1
                                newCellPolygon[k,2] = cellPolygons[i][k][2] - N_y
                            else
                                flag2=2
                                newCellPolygon[k,2] = cellPolygons[i][k][2] + N_y
                            end
                        else
                            newCellPolygon[k,2] = cellPolygons[i][k][2]
                        end

                    end

                    oppositePolygon1[:,1] = newCellPolygon[:,1] .+ N_x
                    oppositePolygon1[:,2] = newCellPolygon[:,2] .+ N_y
                    oppositePolygon2[:,1] = newCellPolygon[:,1] .+ N_x
                    oppositePolygon2[:,2] = newCellPolygon[:,2] .- N_y
                    oppositePolygon3[:,1] = newCellPolygon[:,1] .- N_x
                    oppositePolygon3[:,2] = newCellPolygon[:,2] .+ N_y
                    oppositePolygon4[:,1] = newCellPolygon[:,1] .- N_x
                    oppositePolygon4[:,2] = newCellPolygon[:,2] .- N_y
                    oppositePolygon5[:,1] = newCellPolygon[:,1] .+ N_x
                    oppositePolygon5[:,2] = newCellPolygon[:,2] 
                    oppositePolygon6[:,1] = newCellPolygon[:,1] .- N_x
                    oppositePolygon6[:,2] = newCellPolygon[:,2] 
                    oppositePolygon7[:,1] = newCellPolygon[:,1] 
                    oppositePolygon7[:,2] = newCellPolygon[:,2] .+ N_y
                    oppositePolygon8[:,1] = newCellPolygon[:,1]
                    oppositePolygon8[:,2] = newCellPolygon[:,2] .- N_y

                    # Draw a polygon for cell i with colour determined by concentration
                    poly!(ax, oppositePolygon1, color=sol.u[t][i], colorrange=(minimum(u0_S),maximum(u0_S)), colormap=cmap, strokecolor=:black,
                    strokewidth=0.5)
                    poly!(ax, oppositePolygon2, color=sol.u[t][i], colorrange=(minimum(u0_S),maximum(u0_S)), colormap=cmap, strokecolor=:black,
                    strokewidth=0.5)
                    poly!(ax, oppositePolygon3, color=sol.u[t][i], colorrange=(minimum(u0_S),maximum(u0_S)), colormap=cmap, strokecolor=:black,
                    strokewidth=0.5)
                    poly!(ax, oppositePolygon4, color=sol.u[t][i], colorrange=(minimum(u0_S),maximum(u0_S)), colormap=cmap, strokecolor=:black,
                    strokewidth=0.5)
                    poly!(ax, oppositePolygon5, color=sol.u[t][i], colorrange=(minimum(u0_S),maximum(u0_S)), colormap=cmap, strokecolor=:black,
                    strokewidth=0.5)
                    poly!(ax, oppositePolygon6, color=sol.u[t][i], colorrange=(minimum(u0_S),maximum(u0_S)), colormap=cmap, strokecolor=:black,
                    strokewidth=0.5)
                    poly!(ax, oppositePolygon7, color=sol.u[t][i], colorrange=(minimum(u0_S),maximum(u0_S)), colormap=cmap, strokecolor=:black,
                    strokewidth=0.5)
                    poly!(ax, oppositePolygon8, color=sol.u[t][i], colorrange=(minimum(u0_S),maximum(u0_S)), colormap=cmap, strokecolor=:black,
                    strokewidth=0.5)

                    poly!(ax, newCellPolygon, color=sol.u[t][i], colorrange=(minimum(u0_S),maximum(u0_S)), colormap=cmap, strokecolor=:black, 
                    strokewidth=0.5)

                else # the cell isn't on the periodic boundary
                    if i in cellsTypeA
                        poly!(ax, cellPolygons[i], color=RGB(102/255,178/255,255/255), strokecolor=(:black, 1.0), strokewidth=2)
                    else
                        poly!(ax, cellPolygons[i], color=RGB(255/255,178/255,102/255), strokecolor=(:black, 1.0), strokewidth=2)
                    end
                end

                


            else
                if i in cellsTypeA
                    poly!(ax, cellPolygons[i], color=RGB(102/255,178/255,255/255), strokecolor=(:black, 1.0), strokewidth=2)
                else
                    poly!(ax, cellPolygons[i], color=RGB(255/255,178/255,102/255), strokecolor=(:black, 1.0), strokewidth=2)
                end

            end
            
        end
    end

    # Scatter vertices
    if scatterVertices == 1
        scatter!(ax, Point{2,Float64}.(R), color=:green)
        annotations!(ax, string.(collect(1:length(R))), Point{2,Float64}.(R), color=:green)
    end

    # Scatter edge midpoints
    if scatterEdges == 1
        scatter!(ax, Point{2,Float64}.(edgeMidpoints), color=:blue)
        annotations!(ax, string.(collect(1:length(edgeMidpoints))), Point{2,Float64}.(edgeMidpoints), color=:blue)
    end

    # Scatter cell positions
    if scatterCells == 1
        scatter!(ax, Point{2,Float64}.(cellPositions), color=:red)
        annotations!(ax, string.(collect(1:length(cellPositions))), Point{2,Float64}.(cellPositions), color=:red)
    end

    # Plot resultant forces on vertices (excluding external pressure)
    # NB these forces will be those calculated in the previous integration step and thus will not be exactly up to date for the current vertex positions
    if plotForces == 1
        arrows!(ax, Point{2,Float64}.(R), Vec2f.(sum(F, dims=2)), color=:green)
    end

    if plotEdgeMidpointLinks == 1
        for i = 1:nCells
            for j = 1:cellEdgeCount[i]
                lines!(ax,
                    Point{2,Float64}.([edgeMidpoints[cellEdgeOrders[i][j]],(edgeMidpoints[cellEdgeOrders[i][j]] .+ edgeMidpointLinks[i, cellVertexOrders[i][j]])]),
                    linestyle=:dot,
                    color=:black)
            end
        end
    end

    # Set limits
    reset_limits!(ax)

    # Add frame to movie 
    recordframe!(mov)

    return nothing

end

export visualise

end
