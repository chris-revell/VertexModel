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


function visualise(R, t, fig, ax1, ax2, ax3, cbar1, cbar2, mov, params, matrices, plotCells, scatterEdges, scatterVertices, scatterCells, plotXis,plotForces, plotEdgeMidpointLinks,plotStresses,trackInitialRecoil,plotOrientations)

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
        ξs,
        ξsDir,
        ξsDirScaled,
        cellLabels = matrices
    @unpack initialSystem,
        boundaryType,
        nEdges, 
        nVerts, 
        nCells, 
        L_x,
        L_y,
        k_tracked,
        l_AA,
        l_BB,
        l_AB,
        l_AE,
        l_BE = params

    empty!(ax1)
    empty!(ax2)
    empty!(ax3)
    delete!(cbar1)
    delete!(cbar2)
    


    if boundaryType == "free"
        ax1.title = "t = $(@sprintf("%.3f", t))"
        if plotStresses == 1
            ax2.title = "Plot of AᵢP_effᵢ at t = $(@sprintf("%.3f", t))"
            ax3.title = "Plot of Aᵢξᵢ at t = $(@sprintf("%.3f", t))"
        end
    elseif boundaryType == "periodic"
        ax1.title = "t = $(@sprintf("%.3f", t))"
        ax1.limits = ((0,L_x),(0,L_y))
    
        if plotStresses == 1
            ax2.title = "Plot of AᵢP_effᵢ at t = $(@sprintf("%.3f", t))"
            ax2.limits = ((0,L_x),(0,L_y))
        
            ax3.title = "Plot of Aᵢξᵢ at t = $(@sprintf("%.3f", t))"
            ax3.limits = ((0,L_x),(0,L_y))
        end
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
    # cmap2 = cgrad([
    #     RGB(1.0, 1.0, 1.0),   # white
    #     RGB(1.0, 0.0, 1.0)    # magenta
    # ], 256)

    cmap2 = cgrad([
        RGB(1.0, 1.0, 1.0),   # 0%   White
        RGB(1.0, 0.8, 0.9),   # 25%  Light Pink
        RGB(1.0, 0.0, 0.5),   # 50%  Hot Pink
        RGB(0.6, 0.0, 0.8),   # 75%  Purple
        RGB(0.3, 0.0, 0.5)    # 100% Dark Purple
    ])
    max2 = maximum(A_iξs)
    min2 = 0
    clims2 = (min2, max2)


    # Plot cells
    if plotCells == 1
        cellPolygons = makeCellPolygons(R, params, matrices)
        for i = 1:nCells
            

            if boundaryType == "free"
                if cellLabels[i] == 0
                    poly!(ax1, cellPolygons[i], color=RGB(102/255,178/255,255/255), strokecolor=(:black, 1.0), strokewidth=2)
                else
                    poly!(ax1, cellPolygons[i], color=RGB(255/255,178/255,102/255), strokecolor=(:black, 1.0), strokewidth=2)
                end

                if plotStresses == 1
                    # Plot effective pressures in ax2
                    poly!(ax2, cellPolygons[i], color=A_iP_effs[i], colorrange = clims1, colormap = cmap1,  strokecolor=(:black, 1.0), strokewidth=2)

                    # Plot deviatoric stress in ax3
                    poly!(ax3, cellPolygons[i], color=A_iξs[i], colorrange = clims2, colormap = cmap2,  strokecolor=(:black, 1.0), strokewidth=2)

                end

            elseif boundaryType == "periodic"
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
                    if cellLabels[i] == 0
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

                    if plotStresses == 1
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
                    
                    end    
                else # the cell isn't on the periodic boundary
                    if cellLabels[i] == 0
                        poly!(ax1, cellPolygons[i], color=RGB(102/255,178/255,255/255), strokecolor=(:black, 1.0), strokewidth=2)
                    else
                        poly!(ax1, cellPolygons[i], color=RGB(255/255,178/255,102/255), strokecolor=(:black, 1.0), strokewidth=2)
                    end

                    if plotStresses == 1
                        # Plot effective pressures in ax2
                        poly!(ax2, cellPolygons[i], color=A_iP_effs[i], colorrange = clims1, colormap = cmap1,  strokecolor=(:black, 1.0), strokewidth=2)

                        # Plot deviatoric stress in ax3
                        poly!(ax3, cellPolygons[i], color=A_iξs[i], colorrange = clims2, colormap = cmap2,  strokecolor=(:black, 1.0), strokewidth=2)
                    end
                    
                end


            end
         
        end
    end

    

    # Add the colorbar
    if plotStresses == 1
        cbar1 = Colorbar(fig[2,2],colormap = cmap1,colorrange=clims1, label="AᵢP_effᵢ", width=20,height=Relative(0.6))
        cbar2 = Colorbar(fig[2,3],colormap = cmap2,colorrange=clims2, label="Aᵢξᵢ", width=20,height=Relative(0.6))
    end
    # Scatter vertices
    if scatterVertices == 1
        scatter!(ax1, Point{2,Float64}.(R), color=:green)
        annotations!(ax1, string.(collect(1:length(R))), Point{2,Float64}.(R), color=:green)
    end

    # Scatter edge midpoints
    if scatterEdges == 1
        scatter!(ax1, Point{2,Float64}.(edgeMidpoints), color=:blue, markersize=5)
        annotations!(ax1, string.(collect(1:length(edgeMidpoints))), Point{2,Float64}.(edgeMidpoints), color=:blue, fontsize=4)
    end

    l_vec = [l_AA, l_AB, l_BB, l_AE, l_BE]
    
    for j=1:nEdges 
        edget1Threshold = l_vec[matrices.edgeLabels[j]+1]*0.15
        if matrices.edgeLengths[j] < edget1Threshold
            scatter!(ax1, Point{2,Float64}(edgeMidpoints[j]), color=:red)
            annotations!(ax1, [string(j)], [Point{2,Float64}(edgeMidpoints[j])], color=:red)
        end
    end

    # Scatter cell positions
    if scatterCells == 1
        for i=1:nCells
            p = Point{2,Float64}(cellPositions[i]...)

            # choose color
            col =  :red 

            scatter!(ax1, [p], color=col)
            
            annotations!(ax1, string.(collect(1:length(cellPositions))), Point{2,Float64}.(cellPositions))
        end
    end

    # Plot resultant forces on vertices (excluding external pressure)
    # NB these forces will be those calculated in the previous integration step and thus will not be exactly up to date for the current vertex positions
    if plotForces == 1
        arrows!(ax1, Point{2,Float64}.(R), Vec2f.(sum(F, dims=2)), color=:green)
    end


    if plotXis == 1 
        # for i=1:nCells
        #     barLength = ξs[i]
        #     scatter!(ax3,Point{2,Float64}.(cellPositions))
        # end
        arrows!(ax3,Point{2,Float64}.(cellPositions), 0.5*Vec2f.(ξsDirScaled), color=:black,linewidth = 1,arrowsize=0)
        arrows!(ax3,Point{2,Float64}.(cellPositions), -0.5*Vec2f.(ξsDirScaled), color=:black,linewidth = 1,arrowsize=0)
    end

    if plotOrientations ==1
        arrows!(ax1,Point{2,Float64}.(cellPositions), 0.5*Vec2f.(ξsDir), color=:black,linewidth = 1,arrowsize=0)
        arrows!(ax1,Point{2,Float64}.(cellPositions), -0.5*Vec2f.(ξsDir), color=:black,linewidth = 1,arrowsize=0)
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

    if trackInitialRecoil == 1
        for k in k_tracked 
            scatter!(ax1, Point{2,Float64}(R[k]), color=:red)
        end
    end

    # Set limits
    reset_limits!(ax1)
    if plotStresses == 1
        reset_limits!(ax2)
        reset_limits!(ax3)
    end

    # Add frame to movie 
    recordframe!(mov)

    return cbar1, cbar2

end

export visualise

end
