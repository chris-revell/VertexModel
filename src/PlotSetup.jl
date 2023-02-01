#
#  PlotSetup.jl
#  VertexModel
#
#  Created by Christopher Revell on 31/01/2022.
#
#

module PlotSetup

# Julia packages
using LinearAlgebra
using JLD2
using SparseArrays
using StaticArrays
using UnPack
using DrWatson
using Printf
using FromFile
using CairoMakie

# Local modules
@from "CreateRunDirectory.jl" using CreateRunDirectory
@from "Visualise.jl" using Visualise

function plotSetup(params,matrices,subFolder,folderName)

    # Extract some variables from containers for use below    
    @unpack R = matrices

    # Create plot canvas
    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    fig = Figure(resolution=(2000,1000))
    grid = fig[1,1] = GridLayout()
    ax1 = Axis(grid[1,1],aspect=DataAspect())
    ax2 = Axis(grid[1,2],aspect=DataAspect())
    hidedecorations!(ax1)
    hidespines!(ax1)
    hidedecorations!(ax2)
    hidespines!(ax2)
    xlims!(ax1,min(minimum(first.(R)),minimum(last.(R))), max(maximum(first.(R)),maximum(last.(R))))
    ylims!(ax1,min(minimum(first.(R)),minimum(last.(R))), max(maximum(first.(R)),maximum(last.(R))))
    # Create animation object for visualisation
    mov = VideoStream(fig, framerate=5)
    # Visualise initial system
    visualise(0.0,fig,ax1,ax2,mov,params,matrices)
    save(datadir(subFolder,folderName,"frames","frame$(@sprintf("%03d", 0)).png"),fig)
    
    return fig, ax1, ax2, mov
   
end

export plotSetup 

end
