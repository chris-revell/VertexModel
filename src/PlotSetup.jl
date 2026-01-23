#
#  PlotSetup.jl
#  VertexModel
#
#  Created by Christopher Revell on 31/01/2022.
#
#

module PlotSetup

# Julia packages
using UnPack
using FromFile
using CairoMakie

function plotSetup()

    # Create plot canvas
    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    fig = Figure(size=(1200,600))
    grid = fig[1,1] = GridLayout()
    ax1 = Axis(grid[1,1],aspect=DataAspect())
    ax2 = Axis(grid[1,2],aspect=DataAspect())
    hidedecorations!(ax1)
    hidedecorations!(ax2)
    hidespines!(ax1)
    hidespines!(ax2)

    # Create an empty colorbar in fig[1,3]
    dummy_cmap = cgrad([:white, :white], 2)
    dummy_range = (0.0, 1.0)
    cbar = Colorbar(fig[1,3], colormap=dummy_cmap, colorrange=dummy_range, width=20, height=Relative(0.6))
    

    # Create animation object for visualisation
    mov = VideoStream(fig, framerate=5)
    
    return fig, ax1, ax2, cbar,mov
   
end

export plotSetup 

end
