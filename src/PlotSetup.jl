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
# @from "$(projectdir("src","CreateRunDirectory.jl"))" using CreateRunDirectory
# @from "$(projectdir("src","Visualise.jl"))" using Visualise

function plotSetup(R,params,matrices,subFolder,folderName)

    # Create plot canvas
    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    fig = Figure(resolution=(1000,1000))
    grid = fig[1,1] = GridLayout()
    ax1 = Axis(grid[1,1],aspect=DataAspect())
    hidedecorations!(ax1)
    hidespines!(ax1)
    # Create animation object for visualisation
    mov = VideoStream(fig, framerate=5)
    # Visualise initial system
    # visualise(R,0.0,fig,ax1,mov,params,matrices)
    save(datadir(subFolder,folderName,"frames","frame$(@sprintf("%03d", 0)).png"),fig)
    
    return fig, ax1, mov
   
end

export plotSetup 

end
