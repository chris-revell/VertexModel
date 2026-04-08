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


function visualise(R, t, fig, ax1, ax2, ax3, cbar1, cbar2, mov, params, matrices, plotCells, scatterEdges, scatterVertices, scatterCells, plotForces, plotEdgeMidpointLinks)

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
        P_effs,
        ξs = matrices
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
    empty!(ax3)
    delete!(cbar1)
    delete!(cbar2)
    


    if initialSystem == "new"
        ax1.title = "t = $(@sprintf("%.3f", t))"
        ax2.title = "Plot of P_effᵢ at t = $(@sprintf("%.3f", t))"
        ax3.title = "Plot of ξᵢ at t = $(@sprintf("%.3f", t))"
    else
        ax1.title = "t = $(@sprintf("%.3f", t))"
        ax1.limits = ((0,L_x),(0,L_y))
    
        ax2.title = "Plot of P_effᵢ at t = $(@sprintf("%.3f", t))"
        ax2.limits = ((0,L_x),(0,L_y))
    
        ax3.title = "Plot of ξᵢ at t = $(@sprintf("%.3f", t))"
        ax3.limits = ((0,L_x),(0,L_y))
    end
   
    
    # Set colour limits for P_eff
    A_iP_effs = zeros(nCells)
    A_iP_effs .= cellAreas.*P_effs
    # Generate a colour map for effective pressures: 
    cmap1 = cgrad([
        RGB(0.0, 0.0, 1.0),    # blue
        RGB(1.0, 1.0, 1.0),   # white, zero
        RGB(1.0, 0.0, 0.0)   # red
    ], 256)
    # Alternative colour bar - centered at 0: 
    maxabs1 = maximum(abs.(A_iP_effs))
    clims1 = (-maxabs1, maxabs1)

    # Set colour limits for ξ
    A_iξs = zeros(nCells)
    A_iξs .= cellAreas.*ξs
    # Generate a colour map for ξ:
    cmap2 = cgrad([
        RGB(1.0, 1.0, 1.0),   # white
        RGB(1.0, 0.0, 1.0)    # magenta
    ], 256)
    max2 = maximum(A_iξs)
    min2 = 0
    clims2 = (min2, max2)


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

                # Plot effective pressures in ax2
                poly!(ax2, cellPolygons[i], color=A_iP_effs[i], colorrange = clims1, colormap = cmap1,  strokecolor=(:black, 1.0), strokewidth=2)

                # Plot deviatoric stress in ax3
                poly!(ax3, cellPolygons[i], color=A_iξs[i], colorrange = clims2, colormap = cmap2,  strokecolor=(:black, 1.0), strokewidth=2)

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
                    poly!(ax2, oppositePolygon1, color=A_iP_effs[i], colorrange = clims1, colormap = cmap1, strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax2, oppositePolygon2, color=A_iP_effs[i], colorrange = clims1, colormap = cmap1, strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax2, oppositePolygon3, color=A_iP_effs[i], colorrange = clims1, colormap = cmap1,  strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax2, oppositePolygon4, color=A_iP_effs[i], colorrange = clims1, colormap = cmap1,  strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax2, oppositePolygon5, color=A_iP_effs[i], colorrange = clims1, colormap = cmap1,  strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax2, oppositePolygon6, color=A_iP_effs[i], colorrange = clims1, colormap = cmap1,  strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax2, oppositePolygon7, color=A_iP_effs[i], colorrange = clims1, colormap = cmap1,  strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax2, oppositePolygon8, color=A_iP_effs[i], colorrange = clims1, colormap = cmap1,  strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax2, newCellPolygon, color=A_iP_effs[i], colorrange = clims1, colormap = cmap1,  strokecolor=(:black, 1.0), strokewidth=2)
                    

                    # Plot deviatoric stress in ax3
                    poly!(ax3, oppositePolygon1, color=A_iξs[i], colorrange = clims2, colormap = cmap2, strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax3, oppositePolygon2, color=A_iξs[i], colorrange = clims2, colormap = cmap2, strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax3, oppositePolygon3, color=A_iξs[i], colorrange = clims2, colormap = cmap2,  strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax3, oppositePolygon4, color=A_iξs[i], colorrange = clims2, colormap = cmap2,  strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax3, oppositePolygon5, color=A_iξs[i], colorrange = clims2, colormap = cmap2,  strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax3, oppositePolygon6, color=A_iξs[i], colorrange = clims2, colormap = cmap2,  strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax3, oppositePolygon7, color=A_iξs[i], colorrange = clims2, colormap = cmap2,  strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax3, oppositePolygon8, color=A_iξs[i], colorrange = clims2, colormap = cmap2,  strokecolor=(:black, 1.0), strokewidth=2)
                    poly!(ax3, newCellPolygon, color=A_iξs[i], colorrange = clims2, colormap = cmap2,  strokecolor=(:black, 1.0), strokewidth=2)
                    
                else # the cell isn't on the periodic boundary
                    if i in cellsTypeA
                        poly!(ax1, cellPolygons[i], color=RGB(102/255,178/255,255/255), strokecolor=(:black, 1.0), strokewidth=2)
                    else
                        poly!(ax1, cellPolygons[i], color=RGB(255/255,178/255,102/255), strokecolor=(:black, 1.0), strokewidth=2)
                    end


                    # Plot effective pressures in ax2
                    poly!(ax2, cellPolygons[i], color=A_iP_effs[i], colorrange = clims1, colormap = cmap1,  strokecolor=(:black, 1.0), strokewidth=2)

                    # Plot deviatoric stress in ax3
                    poly!(ax3, cellPolygons[i], color=A_iξs[i], colorrange = clims2, colormap = cmap2,  strokecolor=(:black, 1.0), strokewidth=2)

                    
                end


            end
         
        end
    end

    

    # Add the colorbar
    cbar1 = Colorbar(fig[2,2],colormap = cmap1,colorrange=clims1, label="P_effᵢ", width=20,height=Relative(0.6))
    cbar2 = Colorbar(fig[2,3],colormap = cmap2,colorrange=clims2, label="ξᵢ", width=20,height=Relative(0.6))

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

    return cbar1, cbar2

end

export visualise

end
