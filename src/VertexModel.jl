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
using OrdinaryDiffEq
using LinearAlgebra
using JLD2
using SparseArrays
using StaticArrays
using GLMakie
using CairoMakie
using Printf

# Local modules
@from "$(srcdir("CreateRunDirectory.jl"))" using CreateRunDirectory
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("Initialise.jl"))" using Initialise
@from "$(srcdir("SpatialData.jl"))" using SpatialData
@from "$(srcdir("PlotSetup.jl"))" using PlotSetup
@from "$(srcdir("Model.jl"))" using Model
@from "$(srcdir("T1Transitions.jl"))" using T1Transitions
@from "$(srcdir("TopologyChange.jl"))" using TopologyChange
@from "$(srcdir("Division.jl"))" using Division
@from "$(srcdir("SenseCheck.jl"))" using SenseCheck
@from "$(srcdir("EdgeAblation.jl"))" using EdgeAblation

function vertexModel(;
    initialSystem="new",
    nRows=9,
    nCycles=1,
    realCycleTime=86400.0,
    realTimetMax=nCycles*realCycleTime,
    γ=0.2,
    L₀=0.75,
    A₀=1.0,
    viscousTimeScale=1000.0,
    pressureExternal=0.0,
    peripheralTension=0.0,
    t1Threshold = 0.01,
    surfaceRadius = 20.0,
    surfaceReturnAmplitude = 100.0,
    solver=Tsit5(),
    nBlasThreads=1,
    subFolder="",
    outputTotal=100,
    outputToggle=1,
    frameDataToggle=1,
    frameImageToggle=1,
    printToggle=1,
    videoToggle=1,
    plotCells = 1,
    scatterEdges = 0,
    scatterVertices = 0,
    scatterCells = 0,
    plotForces = 0,
    plotEdgeMidpointLinks = 0,
    setRandomSeed = 0,
    abstol = 1e-7, 
    reltol = 1e-4,
) # All arguments are optional and will be instantiated with these default values if not provided at runtime

    BLAS.set_num_threads(nBlasThreads)

    isodd(nRows)&&(nRows>1)  ? nothing : throw("nRows must be an odd number greater than 1.")

    # Set up initial system, packaging parameters and matrices for system into params and matrices containers from VertexModelContainers.jl
    R, params, matrices = initialise(initialSystem, realTimetMax, γ, L₀, A₀, pressureExternal, viscousTimeScale, outputTotal, t1Threshold, realCycleTime, peripheralTension, setRandomSeed, surfaceRadius, surfaceReturnAmplitude; nRows=nRows)

    # Create fun directory, save parameters, and store directory name for later use.
    if outputToggle == 1
        folderName = createRunDirectory(params,subFolder)
        # Create plot object for later use 
        if frameImageToggle==1 || videoToggle==1
            fig, ax, mov = plotSetup()
        end
    end

    # Set up ODE integrator 
    prob = ODEProblem(model!, R, (0.0, Inf), (params, matrices))
    alltStops = collect(0.0:params.outputInterval:params.tMax)
    integrator = init(prob, solver, tstops=alltStops, abstol=abstol, reltol=reltol, save_on=false, save_start=false, save_end=true)
    outputCounter = [1]

    # Iterate until integrator time reaches max system time 
    while integrator.t <= params.tMax && (integrator.sol.retcode == ReturnCode.Default || integrator.sol.retcode == ReturnCode.Success)
        
        # Output data to file 
        if (integrator.t == alltStops[outputCounter[1]] && outputToggle == 1)
            # Update progress on command line 
            printToggle == 1 ? println("$(@sprintf("%.2f", integrator.t))/$(@sprintf("%.2f", params.tMax)), $(outputCounter[1])/$outputTotal") : nothing
            if frameDataToggle == 1
                # In order to label vertex locations as "R" in data output, create a view of (reference to) integrator.u named R 
                R = @view integrator.u[:]
                jldsave(datadir(folderName, "frameData", "systemData$(@sprintf("%03d", outputCounter[1])).jld2"); matrices, params, R)
            end
            if frameImageToggle == 1 || videoToggle == 1
                # Render visualisation of system and add frame to movie
                visualise(integrator.u, integrator.t, fig, ax, mov, params, matrices)
            end
            # Save still image of this time step 
            frameImageToggle == 1 ? save(datadir(folderName, "frameImages", "frameImage$(@sprintf("%03d", outputCounter[1])).png"), fig) : nothing
            outputCounter[1] += 1
        end

        # Step integrator forwards in time to update vertex positions 
        step!(integrator)

        # Update spatial data (edge lengths, cell areas, etc.) following iteration of the integrator
        spatialData!(integrator.u, params, matrices)

        # Check system for T1 transitions 
        if t1Transitions!(integrator, params, matrices) > 0
            u_modified!(integrator, true)
            # senseCheck(matrices.A, matrices.B; marker="T1") # Check for nonzero values in B*A indicating error in incidence matrices           
            topologyChange!(matrices) # Update system matrices after T1 transition
            spatialData!(integrator.u, params, matrices) # Update spatial data after T1 transition  
            # R = @view integrator.u[:]
        end
        if division!(integrator, params, matrices) > 0
            u_modified!(integrator, true)
            # senseCheck(matrices.A, matrices.B; marker="division") # Check for nonzero values in B*A indicating error in incidence matrices          
            topologyChange!(matrices) # Update system matrices after division 
            spatialData!(integrator.u, params, matrices) # Update spatial data after division 
        end
        # Update cell ages with (variable) timestep used in integration step
        matrices.cellTimeToDivide .-= integrator.dt
        matrices.timeSinceT1 .+= integrator.dt
    end

    # If outputToggle==1, save animation object and save final system matrices
    (outputToggle == 1 && videoToggle == 1) ? save(datadir(folderName, "$(splitpath(folderName)[end]).mp4"), mov) : nothing

    # return matrices, integrator
    return nothing
end

function loadData(relativePath; outputNumber=100)
    data = load(projectdir(relativePath, "frameData", "systemData$(@sprintf("%03d", outputNumber)).jld2"))
    return data["R"], data["matrices"], data["params"]
end

# Ensure code is precompiled
@compile_workload begin
    vertexModel(nCycles=0.01, outputToggle=0, frameDataToggle=0, frameImageToggle=0, printToggle=0, videoToggle=0)
end

export vertexModel
export loadData
export visualise3DInstance
export visualiseSet
export ParametersContainer
export MatricesContainer

end
