#
#  VertexModel.jl
#  VertexModel
#
#  Created by Christopher Revell on 31/01/2021.
#
#
# Input parameters:
# initialSystem    (eg. "single")  String specifying initial system state
# realTimetMax     (eg. 86400.0 )  Real time maximum system run time /seconds
# realCycleTime    (eg. 86400.0 )  Cell cycle time in seconds
# γ                (eg. 0.2     )  Parameter in energy relaxation
# L₀               (eg. 0.75    )  Preferred cell perimeter length
# viscousTimeScale (eg. 20.0    )  Relaxation rate, approx from Sarah's data.
# dt               (eg. 0.01    )  Non dimensionalised time step
# A₀               (eg. 1.0     )  Cell preferred area (1.0 by default)
# pressureExternal (eg. 0.2     )  External pressure applied isotropically to system boundary
# outputTotal      (eg. 20      )  Number of data outputs
# t1Threshold      (eg. 0.01    )  Edge length at which a T1 transition is triggered
# outputToggle     (eg. 1       )  Argument controlling whether data are saved from simulation
# plotToggle       (eg. 1       )  Argument controlling whether plots are produced from simulation
# subFolder        (eg. "Test"  )  Name of subfolder within data directory in which to store results

module VertexModel

# Julia packages
using LinearAlgebra
using JLD2
using SparseArrays
using StaticArrays
using UnPack
using DrWatson
using Makie
using CairoMakie
using DelimitedFiles
using Printf
using FromFile

# Local modules
@from "CreateRunDirectory.jl" using CreateRunDirectory
@from "Visualise.jl" using Visualise
@from "Initialise.jl" using Initialise
@from "Iterate.jl" using Iterate
@from "SpatialData.jl" using SpatialData

function vertexModel(initialSystem,realTimetMax,realCycleTime,γ,L₀,A₀,viscousTimeScale,dt,pressureExternal,t1Threshold,outputTotal,outputToggle,plotToggle;subFolder="")

    @quickactivate

    # Set up initial system, packaging parameters and matrices for system into params and matrices containers from VertexModelContainers.jl
    params,matrices = initialise(initialSystem,realTimetMax,γ,L₀,A₀,pressureExternal,dt,viscousTimeScale,outputTotal,t1Threshold,realCycleTime)

    # Extract some variables from containers for use below
    @unpack tMax, outputInterval = params
    @unpack R, tempR, ΔR, cellAges = matrices

    # Set up output if outputToggle argument == 1
    if outputToggle==1
        # Create fun directory, save parameters, and store directory name for later use.
        folderName = createRunDirectory(params,matrices,subFolder)
        jldsave(datadir(subFolder,folderName,"frames","matrices$(@sprintf("%03d", 0)).jld2");matrices)
        jldsave(datadir(subFolder,folderName,"frames","params$(@sprintf("%03d", 0)).jld2");params)
        if plotToggle==1
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
            # Create animation object for visualisation
            mov = VideoStream(fig, framerate=5)
            # Visualise initial system
            visualise(0.0,fig,ax1,ax2,mov,params,matrices)
            save(datadir(subFolder,folderName,"frames","frame$(@sprintf("%03d", 0)).png"),fig)
        end
    end

    t = 0.0001   # Initial time is very small but slightly above 0 to avoid floating point issues with % operator in output interval calculation
    outCount = 0

    while t<tMax

        # 4 step Runge-Kutta integration
        # 1st step of Runge-Kutta
        iterate!(1,params,matrices,t)
        # 2nd step of Runge-Kutta
        iterate!(2,params,matrices,t)
        # 3rd step of Runge-Kutta
        iterate!(3,params,matrices,t)
        # 4th step of Runge-Kutta
        iterate!(4,params,matrices,t)

        # Result of Runge-Kutta steps
        R .+= ΔR
        t += dt
        cellAges .+= dt

        # Visualise system at every output interval
        if t%outputInterval<dt && outputToggle==1
            outCount += 1
            spatialData!(R,params,matrices)
            jldsave(datadir(subFolder,folderName,"frames","matrices$(@sprintf("%03d", outCount)).jld2");matrices)
            jldsave(datadir(subFolder,folderName,"frames","params$(@sprintf("%03d", outCount)).jld2");params)
            if plotToggle==1
                visualise(t,fig,ax1,ax2,mov,params,matrices)
                save(datadir(subFolder,folderName,"frames","frame$(@sprintf("%03d", outCount)).png"),fig)
            end
            println("$(@sprintf("%.2f", t))/$(@sprintf("%.2f", params.tMax))")
        end
    end

    # If outputToggle==1, save animation object and save final system matrices
    if outputToggle==1
        # Update spatial data from final integration step
        spatialData!(R,params,matrices)
        # Store final system characteristic matrices
        jldsave(datadir(subFolder,folderName,"matricesFinal.jld2");matrices)
        jldsave(datadir(subFolder,folderName,"dataFinal.jld2");params)
        # Save animated gif
        plotToggle==1 ? save(datadir(subFolder,folderName,"$(splitpath(folderName)[end]).mp4"),mov) : nothing
    end

end

export vertexModel

end
