#
#  VertexModel.jl
#  VertexModel
#
#  Created by Christopher Revell on 31/01/2021.
#
#

module VertexModel

# Julia packages
using PrecompileTools
using DrWatson
using FromFile
using DifferentialEquations
using LinearAlgebra
using JLD2
using SparseArrays
using StaticArrays
using CairoMakie
using Printf

# Local modules
@from "CreateRunDirectory.jl" using CreateRunDirectory
@from "Visualise.jl" using Visualise
@from "Initialise.jl" using Initialise
@from "SpatialData.jl" using SpatialData
@from "PlotSetup.jl" using PlotSetup
@from "Model.jl" using Model
@from "T1Transitions.jl" using T1Transitions
@from "TopologyChange.jl" using TopologyChange
@from "Division.jl" using Division
@from "SenseCheck.jl" using SenseCheck

function vertexModel(;
    initialSystem="large",
    realTimetMax=1.5*86400.0,
    realCycleTime=86400.0,
    γ=0.2,
    L₀=0.75,
    A₀=1.0,
    viscousTimeScale=20.0,
    pressureExternal=0.1,
    peripheralTension=0.0,
    t1Threshold=0.05,
    outputTotal=100,
    outputToggle=1,
    frameDataToggle=1,
    frameImageToggle=1,
    printToggle=1,
    videoToggle=1,
    subFolder="",
    solver=Tsit5(),
    nBlasThreads=1
) # All arguments are optional and will be instantiated with these default values if not provided at runtime

    BLAS.set_num_threads(nBlasThreads)

    # Set up initial system, packaging parameters and matrices for system into params and matrices containers from VertexModelContainers.jl
    R,params,matrices = initialise(initialSystem,realTimetMax,γ,L₀,A₀,pressureExternal,viscousTimeScale,outputTotal,t1Threshold,realCycleTime,peripheralTension)

    # Set up output if outputToggle argument == 1
    if outputToggle==1
        # Create fun directory, save parameters, and store directory name for later use.
        folderName = createRunDirectory(R,params,matrices,subFolder)
        if frameImageToggle==1 || videoToggle==1
            # Create plot object for later use 
            fig,ax1,mov=plotSetup(R,params,matrices,subFolder,folderName)
        end
    end

    # Set up ODE integrator 
    prob = ODEProblem(model!,R,(0.0,params.tMax),(params,matrices))
    integrator = init(prob,solver,abstol=1e-7,reltol=1e-4) # Adjust tolerances if you notice unbalanced forces in system that should be at equilibrium

    # Iterate until integrator time reaches max system time 
    while integrator.t<params.tMax
        # Update spatial data (edge lengths, cell areas, etc.)
        spatialData!(integrator.u,params,matrices)
        
        # Output data to file 
        if integrator.t%params.outputInterval<integrator.dt
            # Update progress on command line 
            printToggle==1 ? println("$(@sprintf("%.2f", integrator.t))/$(@sprintf("%.2f", params.tMax)), $(integrator.t*100÷params.tMax)/$outputTotal") : nothing 
            if frameDataToggle==1
                # In order to label vertex locations as "R" in data output, create a view of (reference to) integrator.u named R 
                R = @view integrator.u[:]
                jldsave(datadir("sims",subFolder,folderName,"frames","systemData$(@sprintf("%03d", integrator.t*100÷params.tMax)).jld2");matrices,params,R)
            end
            if frameImageToggle==1 || videoToggle==1
                # Render visualisation of system and add frame to movie
                visualise(integrator.u, integrator.t,fig,ax1,mov,params,matrices)                
            end
            # Save still image of this time step 
            frameImageToggle==1 ? save(datadir("sims",subFolder,folderName,"frames","frame$(@sprintf("%03d", integrator.t*100÷params.tMax)).png"),fig) : nothing
        end

        # Check system for T1 transitions 
        if t1Transitions!(integrator.u,params,matrices)>0
            # senseCheck(matrices.A, matrices.B; marker="T1") # Check for nonzero values in B*A indicating error in incidence matrices           
            topologyChange!(matrices) # Update system matrices after T1 transition  
            spatialData!(integrator.u,params,matrices) # Update spatial data after T1 transition  
        end
        if division!(integrator,params,matrices)>0
            # senseCheck(matrices.A, matrices.B; marker="division") # Check for nonzero values in B*A indicating error in incidence matrices          
            topologyChange!(matrices) # Update system matrices after division 
            spatialData!(integrator.u,params,matrices) # Update spatial data after division 
        end
        
        # Step integrator forwards in time to update vertex positions 
        step!(integrator)

        # Update cell ages with (variable) timestep used in integration step
        matrices.cellAges .+= integrator.dt
    end

    # If outputToggle==1, save animation object and save final system matrices
    if outputToggle==1
        # Update spatial data after final integration step
        spatialData!(integrator.u,params,matrices)
        # Save final system state to file 
        jldsave(datadir("sims",subFolder,folderName,"systemData$outputTotal.jld2");matrices,params,R)        
        # Save movie of simulation if videoToggle==1
        videoToggle==1 ? save(datadir("sims",subFolder,folderName,"$(splitpath(folderName)[end]).mp4"),mov) : nothing
    end

end

# Ensure code is precompiled
@compile_workload begin
    vertexModel(realTimetMax=1000.0,outputToggle=0,frameDataToggle=0,frameImageToggle=0,printToggle=0,videoToggle=0)
end

export vertexModel

end
