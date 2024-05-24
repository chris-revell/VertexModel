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
using DiffEqCallbacks
using DelimitedFiles
using LinearAlgebra
using JLD2
using SparseArrays
using StaticArrays
using CairoMakie
using Printf
using CSV

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
@from "EdgeAblation.jl" using EdgeAblation
@from "Energy.jl" using Energy
@from "Laplacians.jl" using Laplacians

function conditionSteadyState(u, t, integrator)
    # if maximum(norm.(get_du(integrator))) < 1e-3
    #     @show maximum(norm.(get_du(integrator)))
    #     return true
    # else 
    #     @show maximum(norm.(get_du(integrator)))
    #     return false
    # end
    # @show maximum(norm.(get_du(integrator)))
    
    # get_du calculates all gradients as SVectors; 
    # norm.() calculates magnitudes of all gradients as Floats; 
    # maximum() finds biggest gradient
    # Return true if biggest gradient is below threshold 
    #@show maximum(norm.(get_du(integrator)))
    maximum(norm.(get_du(integrator))) < 1e-8   ? true : false
    # Use integrator.opts.abstol as threshold?
end

function affectTerminate!(integrator)
    # if conditionSteadyState() returns true, terminate integrator and pass successful return code
    println("Terminate at steady state")
    terminate!(integrator, ReturnCode.Success)    
end
# Create callback using two user-defined functions above
cb = DiscreteCallback(conditionSteadyState, affectTerminate!)


function vertexModel(;
    initialSystem="large",
    nCycles=1,
    realCycleTime=86400.0,
    realTimetMax=nCycles*realCycleTime,
    γ=0.2,
    L₀=0.75,
    A₀=1.0,
    viscousTimeScale=1000.0,
    pressureExternal=0.0,
    peripheralTension=0.0,
    t1Threshold=0.05,    
    #solver= Vern7(lazy=false),
    solver= DP8(),
    nBlasThreads=1,
    subFolder="",
    outputTotal=100,
    outputToggle=1,
    frameDataToggle=1,
    frameImageToggle=1,
    printToggle=1,
    videoToggle=0,    
    plotCells = 1,
    scatterEdges = 0,
    scatterVertices = 0,
    scatterCells = 0,
    plotForces = 0,
    plotEdgeMidpointLinks = 0,
    setRandomSeed = 0,
) # All arguments are optional and will be instantiated with these default values if not provided at runtime
   
    BLAS.set_num_threads(nBlasThreads)

    # realTimetMax=nCycles*realCycleTime
    # Set up initial system, packaging parameters and matrices for system into params and matrices containers from VertexModelContainers.jl
    R, params, matrices = initialise(initialSystem, realTimetMax, γ, L₀, A₀, pressureExternal, viscousTimeScale, outputTotal, t1Threshold, realCycleTime, peripheralTension, setRandomSeed)

    # Set up output if outputToggle argument == 1
    if outputToggle==1
        # Create fun directory, save parameters, and store directory name for later use.
        folderName = createRunDirectory(R, params, matrices, subFolder)
        if frameImageToggle==1 || videoToggle==1
            # Create plot object for later use 
            fig, ax1, mov = plotSetup(R, params, matrices, subFolder, folderName)
        end
    end


    prob=ODEProblem(model!,R,(0.0,Inf),(params,matrices))

    integrator = init(prob,solver,abstol=1e-10, reltol=1e-8, callback=cb, save_on=false, save_start=false, save_end=true) # Adjust tolerances if you notice unbalanced forces in system that should be at equilibrium

    # Iterate until integrator time reaches max system time 
    while integrator.t < params.tMax && (integrator.sol.retcode == ReturnCode.Default || integrator.sol.retcode == ReturnCode.Success)
        
        # Output data to file 
        if integrator.t % params.outputInterval < integrator.dt
            # Update progress on command line 
            printToggle == 1 ? println("$(@sprintf("%.2f", integrator.t))/$(@sprintf("%.2f", params.tMax)), $(Int64(integrator.t*outputTotal÷params.tMax))/$outputTotal") : nothing
            if frameDataToggle == 1
                # In order to label vertex locations as "R" in data output, create a view of (reference to) integrator.u named R 

                R = @view integrator.u[:]

                jldsave(datadir("sims",subFolder,folderName,"frameData","systemData$(@sprintf("%03d", integrator.t*outputTotal÷params.tMax)).jld2");matrices,params,R)            

                e_tot=energy(params,matrices)
                csvfile=open(datadir("sims",subFolder,folderName,"EnergyTotal.csv"), "a")
                println(csvfile, string(integrator.t*outputTotal÷params.tMax), ",",e_tot)
                close(csvfile)
                print(e_tot ,'\n')
            end
            if frameImageToggle == 1 || videoToggle == 1
                # Render visualisation of system and add frame to movie
                visualise(integrator.u, integrator.t, fig, ax1, mov, params, matrices, plotCells, scatterEdges, scatterVertices, scatterCells, plotForces, plotEdgeMidpointLinks)
            end
            # Save still image of this time step 
            frameImageToggle == 1 ? save(datadir("sims", subFolder, folderName, "frameImages", "frameImage$(@sprintf("%03d", integrator.t*outputTotal÷params.tMax)).png"), fig) : nothing
        end
        # Step integrator forwards in time to update vertex positions 
        step!(integrator)

        # Update spatial data (edge lengths, cell areas, etc.) following iteration of the integrator
        spatialData!(integrator.u, params, matrices)

        # Check system for T1 transitions 
        if t1Transitions!(integrator.u, params, matrices) > 0
            u_modified!(integrator, true)
            # senseCheck(matrices.A, matrices.B; marker="T1") # Check for nonzero values in B*A indicating error in incidence matrices           
            topologyChange!(matrices) # Update system matrices after T1 transition
            spatialData!(integrator.u, params, matrices) # Update spatial data after T1 transition  
        end
       
        if params.nCells < 100
            if division!(integrator,params,matrices)>0
                u_modified!(integrator,true)
                # senseCheck(matrices.A, matrices.B; marker="division") # Check for nonzero values in B*A indicating error in incidence matrices          
                topologyChange!(matrices) # Update system matrices after division 
                spatialData!(integrator.u,params,matrices) # Update spatial data after division 
            end
        end
        #@show params.nCells


       
        
        # Update cell ages with (variable) timestep used in integration step
        matrices.cellTimeToDivide .-= integrator.dt
        matrices.timeSinceT1 .+= integrator.dt
    end

    # If outputToggle==1, save animation object and save final system matrices
    if outputToggle == 1
        # Update spatial data after final integration step
        spatialData!(integrator.u, params, matrices)
        printToggle == 1 ? println("$(@sprintf("%.3f", integrator.t))/$(@sprintf("%.3f", params.tMax)), $(integrator.t*outputTotal÷params.tMax)/$outputTotal") : nothing
        # Save final data file regardless of whether other timepoint data files are saved
        # In order to label vertex locations as "R" in data output, create a view of (reference to) integrator.u named R 
        R = @view integrator.u[:]
        fname=@savename L₀ γ
        jldsave(datadir("sims",subFolder,folderName,"frameData","systemData$(@sprintf("%03d", integrator.t*outputTotal÷params.tMax)).jld2");matrices,params,R)
        jldsave(datadir("sims",subFolder,folderName,"systemDataFinal_$(fname).jld2");matrices,params,R)

        if frameImageToggle==1 || videoToggle==1
            # Render visualisation of system and add frame to movie
            visualise(integrator.u, integrator.t, fig, ax1, mov, params, matrices, plotCells, scatterEdges, scatterVertices, scatterCells, plotForces, plotEdgeMidpointLinks)
        end
        # Save still image of this time step 
        frameImageToggle == 1 ? save(datadir("sims", subFolder, folderName, "frameImages", "frameImage$(@sprintf("%03d", integrator.t*outputTotal÷params.tMax)).png"), fig) : nothing
        # Save movie of simulation if videoToggle==1
        videoToggle == 1 ? save(datadir("sims", subFolder, folderName, "$(splitpath(folderName)[end]).mp4"), mov) : nothing
    end

end

# Ensure code is precompiled
@compile_workload begin
    vertexModel(nCycles=0.01, outputToggle=0, frameDataToggle=0, frameImageToggle=0, printToggle=0, videoToggle=0)
end

export vertexModel

end
