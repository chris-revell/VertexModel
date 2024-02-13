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
# @from "CreateRunDirectory.jl" using CreateRunDirectory
# @from "Visualise.jl" using Visualise

function plotSetup(R,params,matrices,subFolder,folderName)

    # Create plot canvas
    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    fig = Figure(size=(1000,1000))
    grid = fig[1,1] = GridLayout()
    ax1 = Axis(grid[1,1],aspect=DataAspect())
    #hidedecorations!(ax1)
    #hidespines!(ax1)
    xlims!(ax1,-1.25*maximum(abs.([R[x][1] for x in 1:length(R)])),1.25*maximum(abs.([R[x][1] for x in 1:length(R)])))
    ylims!(ax1,-1.25*maximum(abs.([R[x][2] for x in 1:length(R)])),1.25*maximum(abs.([R[x][2] for x in 1:length(R)])))#link to stretch length
    # Create animation object for visualisation
    mov = VideoStream(fig, framerate=5)
    
    return fig, ax1, mov
   
end

export plotSetup 

end
