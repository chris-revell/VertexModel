#
#  VertexModelScript.jl
#  VertexModel
#
#  This script will run all of the components of vertexModel() outside the scope of a function, which can be useful for debugging

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
using Printf
using Dates
using Random
using Distributions
using CircularArrays
using GeometryBasics

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
@from "$(srcdir("VertexModelContainers.jl"))" using VertexModelContainers
@from "$(srcdir("LargeInitialSystem.jl"))" using LargeInitialSystem
@from "$(srcdir("RotationMatrix.jl"))" using RotationMatrix
@from "$(srcdir("ResizeMatrices.jl"))" using ResizeMatrices


initialSystem="large"
nCycles=0.1
realCycleTime=300.0
realTimetMax=nCycles*realCycleTime
γ=0.2
L₀=0.75
A₀=1.0
viscousTimeScale=1000.0
pressureExternal=0.0
peripheralTension=0.0
t1Threshold=0.05
solver=Tsit5()
nBlasThreads=1
subFolder=""
outputTotal=100
outputToggle=1
frameDataToggle=1
frameImageToggle=1
printToggle=1
videoToggle=1
plotCells = 1
scatterEdges = 0
scatterVertices = 0
scatterCells = 0
plotForces = 0
plotEdgeMidpointLinks = 0
setRandomSeed = 0
surfaceRadius = 10

# BLAS.set_num_threads(nBlasThreads)

# realTimetMax=nCycles*realCycleTime
# Set up initial system, packaging parameters and matrices for system into params and matrices containers from VertexModelContainers.jl
# Calculate derived parameters
tMax = realTimetMax / viscousTimeScale  # Non dimensionalised maximum system run time
outputInterval = tMax / outputTotal     # Time interval for storing system data (non dimensionalised)
λ = -2.0 * L₀ * γ
nonDimCycleTime = realCycleTime / viscousTimeScale # Non dimensionalised cell cycle time

# Set random seed value and allocate random number generator
# Random seed set from current unix time, 
# unless non zero value of setRandomSeed is passed, in which case random seed is passed value of setRandomSeed
seed = (setRandomSeed == 0 ? floor(Int64, datetime2unix(now())) : setRandomSeed)
Random.seed!(seed)

R, params, matrices = initialise(initialSystem, realTimetMax, γ, L₀, A₀, pressureExternal, viscousTimeScale, outputTotal, t1Threshold, realCycleTime, peripheralTension, setRandomSeed, surfaceRadius)

# Initial evaluation of matrices based on system topology
topologyChange!(matrices)
spatialData!(R, params, matrices)

# Set up output if outputToggle argument == 1
if outputToggle==1
    # Create fun directory, save parameters, and store directory name for later use.
    folderName = createRunDirectory(R,params,matrices,subFolder)
    if frameImageToggle==1 || videoToggle==1
        # Create plot object for later use 
        fig, ax, mov = plotSetup()
    end
end

# Set up ODE integrator 
prob = ODEProblem(model!, R, (0.0, Inf), (params, matrices))
integrator = init(prob, solver, abstol=1e-7, reltol=1e-4, save_on=false, save_start=false, save_end=true) # Adjust tolerances if you notice unbalanced forces in system that should be at equilibrium

counter=[0]
visualise(integrator.u, integrator.t, fig, ax, mov, params, matrices)
frameImageToggle == 1 ? save(datadir("sims", subFolder, folderName, "frameImages", "frameImage$(@sprintf("%03d", counter[1])).png"), fig) : nothing

#%% 
# Iterate until integrator time reaches max system time 
while integrator.t < params.tMax && (integrator.sol.retcode == ReturnCode.Default || integrator.sol.retcode == ReturnCode.Success)

    # Output data to file 
    if integrator.t % params.outputInterval < integrator.dt
        # Update progress on command line 
        printToggle == 1 ? println("$(@sprintf("%.2f", integrator.t))/$(@sprintf("%.2f", params.tMax)), $(Int64(integrator.t*outputTotal÷params.tMax))/$outputTotal") : nothing
        if frameDataToggle == 1
            # In order to label vertex locations as "R" in data output, create a view of (reference to) integrator.u named R 
            R = @view integrator.u[:]
            jldsave(datadir("sims", subFolder, folderName, "frameData", "systemData$(@sprintf("%03d", integrator.t*outputTotal÷params.tMax)).jld2"); matrices, params, R)
        end
        if frameImageToggle == 1 || videoToggle == 1
            # Render visualisation of system and add frame to movie
            visualise(integrator.u, integrator.t, fig, ax, mov, params, matrices)
        end
        # Save still image of this time step 
        frameImageToggle == 1 ? counter[1]+=1 : nothing 
        frameImageToggle == 1 ? save(datadir("sims", subFolder, folderName, "frameImages", "frameImage$(@sprintf("%03d", counter[1])).png"), fig) : nothing
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
if outputToggle == 1
    # Update spatial data after final integration step
    spatialData!(integrator.u, params, matrices)
    printToggle == 1 ? println("$(@sprintf("%.3f", integrator.t))/$(@sprintf("%.3f", params.tMax)), $(integrator.t*outputTotal÷params.tMax)/$outputTotal") : nothing
    # Save final data file regardless of whether other timepoint data files are saved
    # In order to label vertex locations as "R" in data output, create a view of (reference to) integrator.u named R 
    R = @view integrator.u[:]
    jldsave(datadir("sims", subFolder, folderName, "frameData", "systemData$(@sprintf("%03d", integrator.t*outputTotal÷params.tMax)).jld2"); matrices, params, R)
    if frameImageToggle == 1 || videoToggle == 1
        # Render visualisation of system and add frame to movie
        visualise(integrator.u, integrator.t, fig, ax, mov, params, matrices)
    end
    # Save still image of this time step 
    frameImageToggle == 1 ? save(datadir("sims", subFolder, folderName, "frameImages", "frameImage$(@sprintf("%03d", integrator.t*outputTotal÷params.tMax)).png"), fig) : nothing
    # Save movie of simulation if videoToggle==1
    videoToggle == 1 ? save(datadir("sims", subFolder, folderName, "$(splitpath(folderName)[end]).mp4"), mov) : nothing
end
