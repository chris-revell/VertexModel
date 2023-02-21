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
using DifferentialEquations
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
@from "$(projectdir("src","CreateRunDirectory.jl"))" using CreateRunDirectory
@from "$(projectdir("src","Visualise.jl"))" using Visualise
@from "$(projectdir("src","Initialise.jl"))" using Initialise
@from "$(projectdir("src","Iterate.jl"))" using Iterate
@from "$(projectdir("src","SpatialData.jl"))" using SpatialData
@from "$(projectdir("src","PlotSetup.jl"))" using PlotSetup
@from "$(projectdir("src","Model.jl"))" using Model

function vertexModel(initialSystem,realTimetMax,realCycleTime,γ,L₀,A₀,viscousTimeScale,dt,pressureExternal,peripheralTension,t1Threshold,outputTotal,outputToggle,plotToggle;subFolder="")

    # Set up initial system, packaging parameters and matrices for system into params and matrices containers from VertexModelContainers.jl
    params,matrices = initialise(initialSystem,realTimetMax,γ,L₀,A₀,pressureExternal,dt,viscousTimeScale,outputTotal,t1Threshold,realCycleTime,peripheralTension)

    # Extract some variables from containers for use below
    @unpack tMax, outputInterval = params
    @unpack R, ΔR, cellAges = matrices

    # Set up output if outputToggle argument == 1
    if outputToggle==1
        # Create fun directory, save parameters, and store directory name for later use.
        folderName = createRunDirectory(params,matrices,subFolder)
        jldsave(datadir(subFolder,folderName,"frames","matrices$(@sprintf("%03d", 0)).jld2");matrices)
        jldsave(datadir(subFolder,folderName,"frames","params$(@sprintf("%03d", 0)).jld2");params)
        if plotToggle==1
            fig,ax1,mov=plotSetup(params,matrices,subFolder,folderName)
        end
    end

    outCount = 0

    prob = ODEProblem(model!,matrices.R,(0.0,tMax),[params,matrices])
    integrator = init(prob,Tsit5())

    while integrator.t<tMax
        # spatialData!(R,params,matrices)
        if integrator.t%outputInterval<integrator.dt && outputToggle==1
            outCount += 1
            jldsave(datadir(subFolder,folderName,"frames","matrices$(@sprintf("%03d", outCount)).jld2");matrices)
            jldsave(datadir(subFolder,folderName,"frames","params$(@sprintf("%03d", outCount)).jld2");params)
            if plotToggle==1
                visualise(integrator.t,fig,ax1,mov,params,matrices)
                save(datadir(subFolder,folderName,"frames","frame$(@sprintf("%03d", outCount)).png"),fig)
            end
            println("$(@sprintf("%.2f", integrator.t))/$(@sprintf("%.2f", params.tMax))")
        end

        step!(integrator)

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
