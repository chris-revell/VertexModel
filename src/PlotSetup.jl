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
using GLMakie

function plotSetup() #R,params,matrices,subFolder,folderName)

    # Create plot canvas
    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    fig = Figure(size=(1000,1000))
    ax = Axis3(fig[1,1])
    # hidedecorations!(ax)
    # hidespines!(ax)
    # Create animation object for visualisation
    mov = VideoStream(fig, framerate=5)
    
    return fig, ax, mov
   
end

export plotSetup 

end
