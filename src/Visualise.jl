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
using GeometryBasics: area

# Local modules
@from "OrderAroundCell.jl" using OrderAroundCell
@from "AnalysisFunctions.jl" using AnalysisFunctions


function visualise(R, t, fig, ax1, ax2, cbar, mov, params, matrices, plotCells, scatterEdges, scatterVertices, scatterCells, plotForces, plotEdgeMidpointLinks)

    @unpack cellEdgeCount,
        cellVertexOrders,
        cellEdgeOrders,
        cellPositions,
        cellAreas,
        edgeMidpoints,
        F,
        edgeMidpointLinks,
        μ,
        boundaryCells,
        P_effs = matrices
    @unpack initialSystem,
        nEdges, 
        nVerts, 
        nCells, 
        cellsTypeA, 
        cellsTypeB,
        L_x,
        L_y = params

    empty!(ax1)
    empty!(ax2)
    delete!(cbar)


    ax1.title = "t = $(@sprintf("%.3f", t))"
    
    ax1.limits = ((0,L_x),(0,L_y))

    ax2.title = "Plot of effective pressures at t = $(@sprintf("%.3f", t))"
    
    ax2.limits = ((0,L_x),(0,L_y))
    

    # Compute mean effective cell pressure: 
    A_iP_effs = zeros(nCells)
    A_iP_effs .= cellAreas.*P_effs
    P_eff = sum(A_iP_effs)
    println("Sum A_iP_effs = ",P_eff)

    # Generate a colour map for effective pressures: 
    cmap = cgrad([
        RGB(1.0, 0.0, 0.0),   # red
        RGB(1.0, 1.0, 1.0),   # white, zero
        RGB(0.0, 0.0, 1.0)    # blue
    ], 256)

    # Generate a vector deviations about the mean 
    # ΔP = A_iP_effs .- P_eff
    ΔP = A_iP_effs
    # maxΔP = maximum(ΔP)
    # minΔP = minimum(ΔP)
    # clims = (minΔP, maxΔP)

    # Alternative colour bar - centered at 0: 
    maxabs = maximum(abs.(ΔP))
    clims = (-maxabs, maxabs)


    
    
    # Plot cells
    if plotCells == 1
        cellPolygons = makeCellPolygons(R, params, matrices)
        for i = 1:nCells
            

            if initialSystem == "new"
                if i in cellsTypeA
                    poly!(ax1, cellPolygons[i], color=RGB(102/255,178/255,255/255), strokecolor=(:black, 1.0), strokewidth=2)
                else
                    poly!(ax1, cellPolygons[i], color=RGB(255/255,178/255,102/255), strokecolor=(:black, 1.0), strokewidth=2)
                end
            else 
                # Check whether it is on the periodic boundary: 
                if boundaryCells[i]==1

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
                    flag1=flag2=flag3=flag4=0

                    for k = 1:num_vertices

                        if norm(cellPolygons[i][k][1]-cellPositions[i][1]) > L_x/2
                            
                            if cellPolygons[i][k][1] > cellPositions[i][1]
                                newCellPolygon[k,1] = cellPolygons[i][k][1] - L_x
                                flag1 = 1
                            else
                                newCellPolygon[k,1] = cellPolygons[i][k][1] + L_x
                                flag1 = 2
                            end
                        else
                            newCellPolygon[k,1] = cellPolygons[i][k][1]
                        end

                        if norm(cellPolygons[i][k][2]-cellPositions[i][2]) > L_y/2
                            
                            if cellPolygons[i][k][2] > cellPositions[i][2]
                                flag2 = 1
                                newCellPolygon[k,2] = cellPolygons[i][k][2] - L_y
                            else
                                flag2=2
                                newCellPolygon[k,2] = cellPolygons[i][k][2] + L_y
                            end
                        else
                            newCellPolygon[k,2] = cellPolygons[i][k][2]
                        end

                    end

                    

                    oppositePolygon1[:,1] = newCellPolygon[:,1] .+ L_x
                    oppositePolygon1[:,2] = newCellPolygon[:,2] .+ L_y
                    oppositePolygon2[:,1] = newCellPolygon[:,1] .+ L_x
                    oppositePolygon2[:,2] = newCellPolygon[:,2] .- L_y
                    oppositePolygon3[:,1] = newCellPolygon[:,1] .- L_x
                    oppositePolygon3[:,2] = newCellPolygon[:,2] .+ L_y
                    oppositePolygon4[:,1] = newCellPolygon[:,1] .- L_x
                    oppositePolygon4[:,2] = newCellPolygon[:,2] .- L_y
                    oppositePolygon5[:,1] = newCellPolygon[:,1] .+ L_x
                    oppositePolygon5[:,2] = newCellPolygon[:,2] 
                    oppositePolygon6[:,1] = newCellPolygon[:,1] .- L_x
                    oppositePolygon6[:,2] = newCellPolygon[:,2] 
                    oppositePolygon7[:,1] = newCellPolygon[:,1] 
                    oppositePolygon7[:,2] = newCellPolygon[:,2] .+ L_y
                    oppositePolygon8[:,1] = newCellPolygon[:,1]
                    oppositePolygon8[:,2] = newCellPolygon[:,2] .- L_y

                    # Draw a polygon for cell i with colour determined by cell type
                    if i in cellsTypeA
                        poly!(ax1, oppositePolygon1, color=RGB(102/255,178/255,255/255), strokecolor=(:black, 1.0), strokewidth=2)
                        poly!(ax1, oppositePolygon2, color=RGB(102/255,178/255,255/255), strokecolor=(:black, 1.0), strokewidth=2)
                        poly!(ax1, oppositePolygon3, color=RGB(102/255,178/255,255/255), strokecolor=(:black, 1.0), strokewidth=2)
                        poly!(ax1, oppositePolygon4, color=RGB(102/255,178/255,255/255), strokecolor=(:black, 1.0), strokewidth=2)
                        poly!(ax1, oppositePolygon5, color=RGB(102/255,178/255,255/255), strokecolor=(:black, 1.0), strokewidth=2)
                        poly!(ax1, oppositePolygon6, color=RGB(102/255,178/255,255/255), strokecolor=(:black, 1.0), strokewidth=2)
                        poly!(ax1, oppositePolygon7, color=RGB(102/255,178/255,255/255), strokecolor=(:black, 1.0), strokewidth=2)
                        poly!(ax1, oppositePolygon8, color=RGB(102/255,178/255,255/255), strokecolor=(:black, 1.0), strokewidth=2)
                        poly!(ax1, newCellPolygon, color=RGB(102/255,178/255,255/255), strokecolor=(:black, 1.0), strokewidth=2)
                    else
                        poly!(ax1, oppositePolygon1, color=RGB(255/255,178/255,102/255), strokecolor=(:black, 1.0), strokewidth=2)
                        poly!(ax1, oppositePolygon2, color=RGB(255/255,178/255,102/255), strokecolor=(:black, 1.0), strokewidth=2)
                        poly!(ax1, oppositePolygon3, color=RGB(255/255,178/255,102/255), strokecolor=(:black, 1.0), strokewidth=2)
                        poly!(ax1, oppositePolygon4, color=RGB(255/255,178/255,102/255), strokecolor=(:black, 1.0), strokewidth=2)
                        poly!(ax1, oppositePolygon5, color=RGB(255/255,178/255,102/255), strokecolor=(:black, 1.0), strokewidth=2)
                        poly!(ax1, oppositePolygon6, color=RGB(255/255,178/255,102/255), strokecolor=(:black, 1.0), strokewidth=2)
                        poly!(ax1, oppositePolygon7, color=RGB(255/255,178/255,102/255), strokecolor=(:black, 1.0), strokewidth=2)
                        poly!(ax1, oppositePolygon8, color=RGB(255/255,178/255,102/255), strokecolor=(:black, 1.0), strokewidth=2)
                        poly!(ax1, newCellPolygon, color=RGB(255/255,178/255,102/255), strokecolor=(:black, 1.0), strokewidth=2)
                   
                    end

                    # Plot effective pressures in ax2
                    poly!(ax2, oppositePolygon1, color=ΔP[i], colorrange = clims, colormap = cmap, strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax2, oppositePolygon2, color=ΔP[i], colorrange = clims, colormap = cmap, strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax2, oppositePolygon3, color=ΔP[i], colorrange = clims, colormap = cmap,  strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax2, oppositePolygon4, color=ΔP[i], colorrange = clims, colormap = cmap,  strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax2, oppositePolygon5, color=ΔP[i], colorrange = clims, colormap = cmap,  strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax2, oppositePolygon6, color=ΔP[i], colorrange = clims, colormap = cmap,  strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax2, oppositePolygon7, color=ΔP[i], colorrange = clims, colormap = cmap,  strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax2, oppositePolygon8, color=ΔP[i], colorrange = clims, colormap = cmap,  strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax2, newCellPolygon, color=ΔP[i], colorrange = clims, colormap = cmap,  strokecolor=(:black, 1.0), strokewidth=2)
                    

                else # the cell isn't on the periodic boundary
                    if i in cellsTypeA
                        poly!(ax1, cellPolygons[i], color=RGB(102/255,178/255,255/255), strokecolor=(:black, 1.0), strokewidth=2)
                    else
                        poly!(ax1, cellPolygons[i], color=RGB(255/255,178/255,102/255), strokecolor=(:black, 1.0), strokewidth=2)
                    end


                    # Plot effective pressures in ax2
                    poly!(ax2, cellPolygons[i], color=ΔP[i], colorrange = clims, colormap = cmap,  strokecolor=(:black, 1.0), strokewidth=2)


                    
                end


            end
         
        end
    end

    

    # Add the colorbar
    cbar = Colorbar(fig[1,3],colormap = cmap,colorrange=clims, label="ΔP_eff", width=20,height=Relative(0.6))
    # cbar = Colorbar(fig[1,3],colormap = cmap,colorrange=[-0.2,0.2], label="ΔP_eff", width=20,height=Relative(0.6))
    

    # Scatter vertices
    if scatterVertices == 1
        scatter!(ax1, Point{2,Float64}.(R), color=:green)
        annotations!(ax1, string.(collect(1:length(R))), Point{2,Float64}.(R), color=:green)
    end

    # Scatter edge midpoints
    if scatterEdges == 1
        scatter!(ax1, Point{2,Float64}.(edgeMidpoints), color=:blue)
        annotations!(ax1, string.(collect(1:length(edgeMidpoints))), Point{2,Float64}.(edgeMidpoints), color=:blue)
    end

    # Scatter cell positions
    if scatterCells == 1
        for i=1:nCells
            p = Point{2,Float64}(cellPositions[i]...)

            # choose color
            col = boundaryCells[i] == 1 ? :red : :blue

            scatter!(ax1, [p], color=col)
            
            annotations!(ax1, string.(collect(1:length(cellPositions))), Point{2,Float64}.(cellPositions))
        end
    end

    # Plot resultant forces on vertices (excluding external pressure)
    # NB these forces will be those calculated in the previous integration step and thus will not be exactly up to date for the current vertex positions
    if plotForces == 1
        arrows!(ax1, Point{2,Float64}.(R), Vec2f.(sum(F, dims=2)), color=:green)
    end

    if plotEdgeMidpointLinks == 1
        for i = 1:nCells
            for j = 1:cellEdgeCount[i]
                lines!(ax1,
                    Point{2,Float64}.([edgeMidpoints[cellEdgeOrders[i][j]],(edgeMidpoints[cellEdgeOrders[i][j]] .+ edgeMidpointLinks[i, cellVertexOrders[i][j]])]),
                    linestyle=:dot,
                    color=:black)
            end
        end
    end

    # Set limits
    reset_limits!(ax1)

    # Add frame to movie 
    recordframe!(mov)

    return cbar

end

export visualise

end
