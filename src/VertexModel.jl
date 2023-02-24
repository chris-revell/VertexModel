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
# A₀               (eg. 1.0     )  Cell preferred area (1.0 by default)
# pressureExternal (eg. 0.2     )  External pressure applied isotropically to system boundary
# outputTotal      (eg. 20      )  Number of data outputs
# t1Threshold      (eg. 0.01    )  Edge length at which a T1 transition is triggered
# outputToggle     (eg. 1       )  Argument controlling whether data are saved from simulation
# plotToggle       (eg. 1       )  Argument controlling whether plots are produced from simulation
# subFolder        (eg. "Test"  )  Name of subfolder within data directory in which to store results

module VertexModel

# Julia packages
using DrWatson
using FromFile
using DifferentialEquations
using LinearAlgebra
using JLD2
using SparseArrays
using StaticArrays
using CairoMakie
using Printf
# using Makie
# using DelimitedFiles


# Local modules
@from "$(projectdir("src","CreateRunDirectory.jl"))" using CreateRunDirectory
@from "$(projectdir("src","Visualise.jl"))" using Visualise
@from "$(projectdir("src","Initialise.jl"))" using Initialise
@from "$(projectdir("src","SpatialData.jl"))" using SpatialData
@from "$(projectdir("src","PlotSetup.jl"))" using PlotSetup
@from "$(projectdir("src","Model.jl"))" using Model
@from "$(projectdir("src","T1Transitions.jl"))" using T1Transitions
@from "$(projectdir("src","TopologyChange.jl"))" using TopologyChange
@from "$(projectdir("src","Division.jl"))" using Division
@from "$(projectdir("src","SenseCheck.jl"))" using SenseCheck

function vertexModel(initialSystem,realTimetMax,realCycleTime,γ,L₀,A₀,viscousTimeScale,pressureExternal,peripheralTension,t1Threshold,outputTotal,outputToggle,plotToggle;subFolder="")

    # Set up initial system, packaging parameters and matrices for system into params and matrices containers from VertexModelContainers.jl
    R,params,matrices = initialise(initialSystem,realTimetMax,γ,L₀,A₀,pressureExternal,viscousTimeScale,outputTotal,t1Threshold,realCycleTime,peripheralTension)

    # Set up output if outputToggle argument == 1
    if outputToggle==1
        # Create fun directory, save parameters, and store directory name for later use.
        folderName = createRunDirectory(R,params,matrices,subFolder)
        if plotToggle==1
            fig,ax1,mov=plotSetup(R,params,matrices,subFolder,folderName)
        end
    end

    # Set up ODE integrator 
    prob = ODEProblem(model!,R,(0.0,params.tMax),[params,matrices])
    integrator = init(prob,Tsit5())

    while integrator.t<params.tMax
        spatialData!(integrator.u,params,matrices)
        
        if integrator.t%params.outputInterval<integrator.dt && outputToggle==1
            jldsave(datadir(subFolder,folderName,"frames","matrices$(@sprintf("%03d", integrator.t*100÷params.tMax)).jld2");matrices)
            jldsave(datadir(subFolder,folderName,"frames","params$(@sprintf("%03d", integrator.t*100÷params.tMax)).jld2");params)
            if plotToggle==1
                visualise(integrator.u, integrator.t,fig,ax1,mov,params,matrices)
                save(datadir(subFolder,folderName,"frames","frame$(@sprintf("%03d", integrator.t*100÷params.tMax)).png"),fig)
            end
            println("$(@sprintf("%.2f", integrator.t))/$(@sprintf("%.2f", params.tMax)), $(integrator.t*100÷params.tMax)/$outputTotal")
        end

        if t1Transitions!(integrator.u,params,matrices)>0
            senseCheck(matrices.A, matrices.B; marker="T1")
            topologyChange!(matrices)
            spatialData!(integrator.u,params,matrices)
        end
        if params.nCells <22
            if division!(integrator,params,matrices)>0
                senseCheck(matrices.A, matrices.B; marker="division")
                topologyChange!(matrices)
                spatialData!(integrator.u,params,matrices)
            end
        end
        step!(integrator)

        matrices.cellAges .+= integrator.dt
    end

    # If outputToggle==1, save animation object and save final system matrices
    if outputToggle==1
        # Update spatial data from final integration step
        spatialData!(integrator.u,params,matrices)
        # Store final system characteristic matrices
        jldsave(datadir(subFolder,folderName,"matricesFinal.jld2");matrices)
        jldsave(datadir(subFolder,folderName,"dataFinal.jld2");params)
        # Save animated gif
        plotToggle==1 ? save(datadir(subFolder,folderName,"$(splitpath(folderName)[end]).mp4"),mov) : nothing
    end

end

export vertexModel

end
