#
#  PlotSetup.jl
#  VertexModel
#
#  Created by Christopher Revell on 31/01/2022.
#  Edited by Charlotte Taylor Barca
#

module PlotSetup

# Julia packages
using UnPack
using FromFile
using CairoMakie

function plotSetup()

    # Create plot canvas
    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    fig = Figure(size=(1800,1200))
    grid = fig[1,1] = GridLayout()
    ax1 = Axis(grid[1,1],aspect=DataAspect())
    ax2 = Axis(grid[1,2],aspect=DataAspect())
    ax3 = Axis(grid[1,3],aspect=DataAspect())
    hidedecorations!(ax1)
    hidedecorations!(ax2)
    hidedecorations!(ax3)
    hidespines!(ax1)
    hidespines!(ax2)
    hidespines!(ax3)

    # Create an empty colorbar in fig[1,3]
    dummy_cmap = cgrad([:white, :white], 2)
    dummy_range = (0.0, 1.0)
    cbar1 = Colorbar(fig[2,2], colormap=dummy_cmap, colorrange=dummy_range, width=20, height=Relative(0.6))
    cbar2 = Colorbar(fig[2,3], colormap=dummy_cmap, colorrange=dummy_range, width=20, height=Relative(0.6))

    # Create animation object for visualisation
    mov = VideoStream(fig, framerate=5)
    
    return fig, ax1, ax2, ax3, cbar1, cbar2, mov
   
end

function stressPlotSetup()

    # Initialise figure for 3 subplots 
    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    stressFig = Figure(size=(1800,600))

    # Initialise a figure for tracking sum of P_effsA_i: 
    gridPeff = stressFig[1,1] = GridLayout()
    axPeff = Axis(gridPeff[1,1],aspect=1)
    axPeff.title = "Sum of AᵢP_effᵢ over time"
    axPeff.xlabel = "t"
    axPeff.ylabel = "ΣᵢAᵢP_effᵢ"

    # Initialise a figure for tracking sum of U_i: 
    gridU =  stressFig[1,2] = GridLayout()
    axU = Axis(gridU[1,1],aspect=1)
    axU.title = "Sum of Uᵢ over time"
    axU.xlabel = "t"
    axU.ylabel = "ΣᵢUᵢ"

    # Initialise figure to track sum of A_iξ_i:
    gridξ =  stressFig[1,3] = GridLayout()
    axξ = Axis(gridξ[1,1],aspect=1)
    axξ.title = "Sum of Aᵢξᵢ over time"
    axξ.xlabel = "t"
    axξ.ylabel = "ΣᵢAᵢξᵢ"

    return stressFig, axPeff, axU, axξ

end # end function

function recoilVelocityPlotSetup()
    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    recoilFig = Figure(size=(600,600))

    # Initialise a figure for tracking sum of P_effsA_i: 
    gridRecoil = recoilFig[1,1] = GridLayout()
    axVelocity = Axis(gridRecoil[1,1],aspect=1)
    axVelocity.title = "Initial recoil velocity"
    axVelocity.xlabel = "t"
    axVelocity.ylabel = "|r₁-r₂|"

    return recoilFig, axVelocity
end

export plotSetup, stressPlotSetup, recoilVelocityPlotSetup

end # end module 
