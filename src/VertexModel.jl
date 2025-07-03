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
using CairoMakie
using Printf
using DiffEqCallbacks

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
@from "Callbacks.jl" using Callbacks

function vertexModel(;
    initialSystem = "new",
    nRows = 9,
    nCycles = 1,
    realCycleTime = 86400.0,
    realTimetMax = nCycles*realCycleTime,
    γ = 0.2,
    L₀ = 0.75,
    A₀ = 1.0,
    viscousTimeScale = 1000.0,
    pressureExternal = 0.0,
    peripheralTension = 0.0,
    t1Threshold = 0.05,
    divisionToggle = 1,
    solver = Tsit5(),
    nBlasThreads = 1,
    subFolder = "",
    outputTotal = 100,
    outputToggle = 1,
    frameDataToggle = 1,
    frameImageToggle = 1,
    printToggle = 1,
    videoToggle = 1,
    plotCells = 1,
    scatterEdges = 0,
    scatterVertices = 0,
    scatterCells = 0,
    plotForces = 0,
    plotEdgeMidpointLinks = 0,
    setRandomSeed = 0,
    abstol = 1e-7, 
    reltol = 1e-4,
    energyModel = "log",
    vertexWeighting = 1,
    R_in = spzeros(2),
    A_in = spzeros(2),
    B_in = spzeros(2), 
    termSteadyState = false,
) # All arguments are optional and will be instantiated with these default values if not provided at runtime

    BLAS.set_num_threads(nBlasThreads)

    # Set up initial system, packaging parameters and matrices for system into params and matrices containers from VertexModelContainers.jl
    u0, params, matrices = initialise(initialSystem = initialSystem,
        realTimetMax = realTimetMax,
        γ = γ,
        L₀ = L₀,
        A₀ = A₀,
        pressureExternal = pressureExternal,
        viscousTimeScale = viscousTimeScale,
        outputTotal = outputTotal,
        t1Threshold = t1Threshold,
        realCycleTime = realCycleTime,
        peripheralTension = peripheralTension,
        setRandomSeed = setRandomSeed,
        nRows = nRows,
        energyModel = energyModel,
        vertexWeighting = vertexWeighting,
        R_in = R_in,
        A_in = A_in,
        B_in = B_in,
    )

    # Create directory in which to store date. Save parameters and store directory name for later use.
    if outputToggle == 1
        folderName = createRunDirectory(params,subFolder)
        # Create plot object for later use 
        if frameImageToggle==1 || videoToggle==1
            fig, ax, mov = plotSetup()
        end
    end

    # Set up ODE integrator   
    prob = ODEProblem(model!,
            u0,
            (0.0, (termSteadyState ? Inf : params.tMax)),
            (params, matrices),
        )
    alltStops = collect(0.0:params.outputInterval: (termSteadyState ? 100.0.*params.tMax : params.tMax)) # Time points that the solver will be forced to land at during integration
    integrator = init(prob,
            solver,
            tstops=alltStops,
            abstol=abstol,
            reltol=reltol,
            save_on=false,
            save_start=false,
            save_end=true,
            callback = DiscreteCallback(termSteadyState ? conditionSteadyState : conditiontMax, affectTerminate!), #(termSteadyState ? cbSS : cbtMax),
        )  
    outputCounter = [1]

    # Iterate until integrator terminates according to specified callback 
    while (integrator.sol.retcode == ReturnCode.Default || integrator.sol.retcode == ReturnCode.Success) && integrator.sol.retcode!=:Terminate
        
        # Reinterpret state vector as a vector of SVectors 
        R = reinterpret(SVector{2,Float64}, integrator.u)
        # Note that reinterpreting accesses the same underlying data, 
        # so changes to R will update integrator.u and vice versa 

        # Output data to file 
        if integrator.t == alltStops[outputCounter[1]]
            termSteadyState  ? (@show maximum(norm.(get_du(integrator)))) : nothing
            # Update progress on command line 
            printToggle == 1 ? println("$(@sprintf("%.2f", integrator.t))/$(@sprintf("%.2f", params.tMax)), $(outputCounter[1])/$outputTotal") : nothing            
            if frameDataToggle == 1
                # Save system data to file 
                jldsave(datadir(folderName, "frameData", "systemData$(@sprintf("%03d", outputCounter[1])).jld2"); matrices, params, R)
            end
            if frameImageToggle == 1 || videoToggle == 1
                # Render visualisation of system and add frame to movie
                visualise(R, integrator.t, fig, ax, mov, params, matrices, plotCells, scatterEdges, scatterVertices, scatterCells, plotForces, plotEdgeMidpointLinks)
            end
            # Save still image of this time step 
            frameImageToggle == 1 ? save(datadir(folderName, "frameImages", "frameImage$(@sprintf("%03d", outputCounter[1])).png"), fig) : nothing
            outputCounter[1] += 1
        end

        # Step integrator forwards in time to update vertex positions 
        step!(integrator)

        # Update spatial data (edge lengths, cell areas, etc.) following iteration of the integrator
        spatialData!(R, params, matrices)

        # Check system for T1 transitions 
        if t1Transitions!(integrator, params, matrices) > 0
            u_modified!(integrator, true)
            # senseCheck(matrices.A, matrices.B; marker="T1") # Check for nonzero values in B*A indicating error in incidence matrices           
            topologyChange!(matrices) # Update system matrices after T1 transition
            spatialData!(R, params, matrices) # Update spatial data after T1 transition  
        end
        if divisionToggle==1
            if division!(integrator, params, matrices) > 0
                u_modified!(integrator, true)
                # senseCheck(matrices.A, matrices.B; marker="division") # Check for nonzero values in B*A indicating error in incidence matrices          
                topologyChange!(matrices) # Update system matrices after division 
                spatialData!(R, params, matrices) # Update spatial data after division 
            end
        end
        # Update cell ages with (variable) timestep used in integration step
        matrices.cellTimeToDivide .-= integrator.dt
        matrices.timeSinceT1 .+= integrator.dt
    end

    # If outputToggle==1, save animation object and save final system matrices
    (outputToggle == 1 && videoToggle == 1) ? save(datadir(folderName, "$(splitpath(folderName)[end]).mp4"), mov) : nothing

    return integrator
end

# Function to load previously saved simulation data 
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

end
