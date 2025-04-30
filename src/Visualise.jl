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
@from "PlotSetup.jl" using PlotSetup

function visualise(R, t, fig, ax, mov, params, matrices, plotCells, scatterEdges, scatterVertices, scatterCells, plotForces, plotEdgeMidpointLinks)

    @unpack cellEdgeCount,
        cellVertexOrders,
        cellEdgeOrders,
        cellPositions,
        edgeMidpoints,
        F,
        edgeMidpointLinks,
        μ,
        cellLineage = matrices
    @unpack nEdges,
        nVerts,
        nCells = params

    empty!(ax)

    # ax.title = "t = $(@sprintf("%.3f", t))"

    # # Plot cells
    # if plotCells == 1
    #     cellPolygons = makeCellPolygons(R, params, matrices)
    #     for i = 1:nCells
    #         poly!(ax, cellPolygons[i], color=(getRandomColor(i), 0.5), strokecolor=(:black, 1.0), strokewidth=2)
    #     end
    # end

    ax.title = "t = $(@sprintf("%.2f", t))"

    # Plot cells
    # if plotCells == 1
    #     cellPolygons = makeCellPolygons(R,params,matrices)
    #     for i=1:nCells
    #         #poly!(ax,cellPolygons[i],color=cellAreas[i],colormap=:viridis,colorrange=(minimum(cellAreas)-1e-6, maximum(cellAreas)+1e-6),strokecolor=(:black,1.0),strokewidth=1)
    #         poly!(ax,cellPolygons[i], color=cellEdgeCount[i], colorrange=(3, 10),colormap=cgrad(ColorSchemes.jet, 8, categorical=true),strokecolor=:black, strokewidth=1)

    #     end
    # end
    
    # cbar=Colorbar(fig[1,2],limits=(3,10),colormap=cgrad(ColorSchemes.jet, 8, categorical=true),flipaxis=true)
    # cbar.ticks = ([3+0.5*(7/8), 3+1.5*(7/8),  3+2.5*(7/8), 3+3.5*(7/8),  3+4.5*(7/8),  3+5.5*(7/8),  3+6.5*(7/8),  3+7.5*(7/8)], ["3", "4", "5","6", "7", "8", "9", "10"])

    if plotCells == 1
        cellPolygons = makeCellPolygons(R,params,matrices)
        for i=1:nCells
            #poly!(ax,cellPolygons[i],color=cellAreas[i],colormap=:viridis,colorrange=(minimum(cellAreas)-1e-6, maximum(cellAreas)+1e-6),strokecolor=(:black,1.0),strokewidth=1)
            poly!(ax,cellPolygons[i], color=cellLineage[i], colorrange=(1, maximum(cellLineage)),colormap=cgrad(ColorSchemes.turbo, maximum(cellLineage), categorical=true),strokecolor=:black, strokewidth=1)

        end
    end
    
    cbar=Colorbar(fig[1,2],limits=(1, maximum(cellLineage)),colormap=cgrad(ColorSchemes.turbo, maximum(cellLineage), categorical=true),flipaxis=true)
    #cbar.ticks = ([3+0.5*(7/8), 3+1.5*(7/8),  3+2.5*(7/8), 3+3.5*(7/8),  3+4.5*(7/8),  3+5.5*(7/8),  3+6.5*(7/8),  3+7.5*(7/8)], ["3", "4", "5","6", "7", "8", "9", "10"])

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

function visualise_final_state(R, t,params, matrices, folderName)

    @unpack cellEdgeCount,
        cellVertexOrders,
        cellEdgeOrders,
        cellPositions,
        edgeMidpoints,
        F,
        edgeMidpointLinks,
        μ,
        cellLineage,
        cellGeneration = matrices
    @unpack nEdges,
        nVerts,
        nCells,
        viscousTimeScale,L₀, γ, 
        realCycleTime = params


    fname=@savename L₀ γ realCycleTime
     
    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")

    #Plot cells by edge count

    fig = Figure(size=(1000,1000))
    grid = fig[1,1] = GridLayout()
    ax = Axis(grid[1,1],aspect=DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)

    ax.title = "t = $(@sprintf("%.2f", t*viscousTimeScale))"

    #Plot cells by edge count

    cellPolygons = makeCellPolygons(R,params,matrices)
    for i=1:nCells
        #poly!(ax,cellPolygons[i],color=cellAreas[i],colormap=:viridis,colorrange=(minimum(cellAreas)-1e-6, maximum(cellAreas)+1e-6),strokecolor=(:black,1.0),strokewidth=1)
        poly!(ax,cellPolygons[i], color=cellEdgeCount[i], colorrange=(3, 10),colormap=cgrad(ColorSchemes.jet, 8, categorical=true),strokecolor=:black, strokewidth=1)

    end
    
    cbar=Colorbar(fig[1,2],limits=(3,10),colormap=cgrad(ColorSchemes.jet, 8, categorical=true),flipaxis=true)
    cbar.ticks = ([3+0.5*(7/8), 3+1.5*(7/8),  3+2.5*(7/8), 3+3.5*(7/8),  3+4.5*(7/8),  3+5.5*(7/8),  3+6.5*(7/8),  3+7.5*(7/8)], ["3", "4", "5","6", "7", "8", "9", "10"])
    save(datadir(folderName, "Cell_Edge_count_Final$(fname).png"), fig)
    
    fig = Figure(size=(1000,1000))
    grid = fig[1,1] = GridLayout()
    ax = Axis(grid[1,1],aspect=DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)

    ax.title = "t = $(@sprintf("%.2f", t*viscousTimeScale))"


    cellPolygons = makeCellPolygons(R,params,matrices)
    for i=1:nCells
        #poly!(ax,cellPolygons[i],color=cellAreas[i],colormap=:viridis,colorrange=(minimum(cellAreas)-1e-6, maximum(cellAreas)+1e-6),strokecolor=(:black,1.0),strokewidth=1)
        poly!(ax,cellPolygons[i], color=cellLineage[i], colorrange=(1, maximum(cellLineage)),colormap=cgrad(ColorSchemes.turbo, maximum(cellLineage), categorical=true),strokecolor=:black, strokewidth=1)
    end

    cbar=Colorbar(fig[1,2],limits=(1, maximum(cellLineage)),colormap=cgrad(ColorSchemes.turbo, maximum(cellLineage),categorical=true),flipaxis=true)

    save(datadir(folderName, "Cell_Lineage_Final$(fname).png"), fig)
    
    fig = Figure(size=(1000,1000))
    grid = fig[1,1] = GridLayout()
    ax = Axis(grid[1,1],aspect=DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)

    ax.title = "t = $(@sprintf("%.2f", t*viscousTimeScale))"


    cellPolygons = makeCellPolygons(R,params,matrices)
    for i=1:nCells
        #poly!(ax,cellPolygons[i],color=cellAreas[i],colormap=:viridis,colorrange=(minimum(cellAreas)-1e-6, maximum(cellAreas)+1e-6),strokecolor=(:black,1.0),strokewidth=1)
        poly!(ax,cellPolygons[i], color=cellGeneration[i], colorrange=(minimum(cellGeneration), maximum(cellGeneration)),colormap=cgrad(ColorSchemes.tableau_jewel_bright, maximum(cellGeneration), categorical=true),strokecolor=:black, strokewidth=1)

    end
    cbar=Colorbar(fig[1,2],limits=(minimum(cellGeneration),maximum(cellGeneration)),colormap=cgrad(ColorSchemes.tableau_jewel_bright, maximum(cellGeneration), categorical=true),flipaxis=true)

    save(datadir(folderName, "Cell_Generation_Final$(fname).png"), fig)


    return nothing
end    



export visualise, visualise_final_state

end
