#
#  PlotSetup.jl
#  VertexModel
#
#

module PlotSetup

# Julia packages
using UnPack
using FromFile
# using GLMakie
using CairoMakie    

function plotSetup()

    CairoMakie.activate!()
    # Create plot canvas
    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    fig = Figure(size=(1000,1000))
    ax = Axis3(fig[1,1], aspect=:data)
    # hidedecorations!(ax)
    # hidespines!(ax)

    # Create animation object for visualisation
    mov = VideoStream(fig, framerate=5)
    
    return fig, ax, mov
   
end

export plotSetup 

end
