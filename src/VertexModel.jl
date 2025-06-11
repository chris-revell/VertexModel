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
@from "Stretch.jl" using Stretch

###
#Using vertex weighted log model
# For no division Vern7 works well, with lazy=false to acount for T1Transitions
# TanYam7 works with division and convergest to SS for 100 cells relatively quickly (ncycles=5, 26 time points, reltol=1e-7), also ok for 200 but with lower tolerances (1e-9, 1e-8)
# Dp8 worked well for 50 cells but didn't converge to SS for 100 cells in 5 cycles
# Tsit5 does not reach steady state
#    
###

function vertexModel(;
    initialSystem="hex",
    nRows=11,
    nCycles=1,
    realCycleTime=14400.0, #4hrs
    realTimetMax=nCycles*realCycleTime,
    γ=0.172, #from ANB parameter inference Xenopus animal caps
    L₀=0.753, #from ANB parameter inference Xenopus animal caps
    A₀=1.0,
    viscousTimeScale=1000.0,
    pressureExternal=0.0,
    peripheralTension=0.0,
    t1Threshold=0.01,
    #solver=TanYam7(),
    #solver=Tsit5(),
    solver=Vern7(lazy=false),
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
    abstol = 1e-10, 
    reltol = 1e-8,
    modelChoice="quadratic",
    vertexWeighting=0,
    stretchType="none", 
    realStretchTime=0,
    λs=0,
    κ=1,
    maxCells=1,
) # All arguments are optional and will be instantiated with these default values if not provided at runtime

    BLAS.set_num_threads(nBlasThreads)

    isodd(nRows)&&(nRows>1)  ? nothing : throw("nRows must be an odd number greater than 1.")

    # Set up initial system, packaging parameters and matrices for system into params and matrices containers from VertexModelContainers.jl
    u0, params, matrices = initialise(initialSystem, realTimetMax, γ, L₀, A₀, pressureExternal, viscousTimeScale, outputTotal, t1Threshold, realCycleTime, peripheralTension, setRandomSeed; nRows=nRows,modelChoice=modelChoice,
    vertexWeighting=vertexWeighting, stretchType=stretchType, realStretchTime=realStretchTime, λs=λs, κ=κ )

    # Create directory in which to store date. Save parameters and store directory name for later use.
    folderName = createRunDirectory(params,subFolder)
    fname=@savename L₀ γ λs realStretchTime κ

    R0 = reinterpret(SVector{2,Float64}, u0)
    jldsave(datadir(folderName,"systemDataInitial_$(fname).jld2");matrices,params,R0)

    # if outputToggle == 1
    #     #Create plot object for later use 

    #     if frameImageToggle==1 || videoToggle==1
            fig, ax, mov = plotSetup()
    #     end
    # end
    visualise(R0, 0, fig, ax, mov, params, matrices, plotCells, scatterEdges, scatterVertices, scatterCells, plotForces, plotEdgeMidpointLinks)
    save(datadir(folderName, "systemDataInitial_$(fname).png"), fig)

    # Set up ODE integrator 

    ###Need to set up tstops with max stretch time and then regular intervals relative to stretch

    prob = ODEProblem(model!, u0, (0.0, Inf), (params, matrices))
    alltStops = collect(0.0:params.outputInterval:params.tMax) # Time points that the solver will be forced to land at during integration
    push!(alltStops,params.tStretch )
    integrator = init(prob, solver, tstops=alltStops, abstol=abstol, reltol=reltol, save_on=false, save_start=false, save_end=true, callback=TerminateSteadyState(min_t=params.tStretch+1))
    outputCounter = [1]

    # Iterate until integrator time reaches max system time 
    while integrator.t <= params.tMax && (integrator.sol.retcode == ReturnCode.Default || integrator.sol.retcode == ReturnCode.Success)
        
        # Reinterpret state vector as a vector of SVectors 
        R = reinterpret(SVector{2,Float64}, integrator.u)
        # Note that reinterpreting accesses the same underlying data, so changes to R will update integrator.u and vice versa 

        # Output data to file 
        if integrator.t == alltStops[outputCounter[1]]
            # Update progress on command line 
            printToggle == 1 ? println("$(@sprintf("%.2f", integrator.t))/$(@sprintf("%.2f", params.tMax)), $(outputCounter[1])/$outputTotal") : nothing            
            if frameDataToggle == 1
                # Save system data to file 
                jldsave(datadir(folderName, "frameData", "systemData$(@sprintf("%03d", outputCounter[1]-1)).jld2"); matrices, params, R)
            end
            if frameImageToggle == 1 || videoToggle == 1
                # Render visualisation of system and add frame to movie
                visualise(R, integrator.t, fig, ax, mov, params, matrices, plotCells, scatterEdges, scatterVertices, scatterCells, plotForces, plotEdgeMidpointLinks)
            end
            # Save still image of this time step 
            frameImageToggle == 1 ? save(datadir(folderName, "frameImages", "frameImage$(@sprintf("%03d", outputCounter[1]-1)).png"), fig) : nothing
            outputCounter[1] += 1
        end
        
        if integrator.t==params.tStretch
            jldsave(datadir(folderName,"systemDataFullStretch_$(fname).jld2");matrices,params,R)

            visualise(R, integrator.t, fig, ax, mov, params, matrices, plotCells, scatterEdges, scatterVertices, scatterCells, plotForces, plotEdgeMidpointLinks)
            save(datadir(folderName, "FullStretch_$(fname).png"), fig)
        end

        # Step integrator forwards in time to update vertex positions 
        step!(integrator)


        # Update spatial data (edge lengths, cell areas, etc.) following iteration of the integrator
        spatialData!(R, params, matrices)
        @show minimum(matrices.edgeLengths)
        # Check system for T1 transitions 
        if t1Transitions!(integrator, params, matrices) > 0
            u_modified!(integrator, true)
            # senseCheck(matrices.A, matrices.B; marker="T1") # Check for nonzero values in B*A indicating error in incidence matrices           
            topologyChange!(matrices) # Update system matrices after T1 transition
            spatialData!(R, params, matrices) # Update spatial data after T1 transition  
        end
        if params.nCells < maxCells
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
    R = reinterpret(SVector{2,Float64}, integrator.u)
    jldsave(datadir(folderName,"systemDataFinal_$(fname).jld2");matrices,params,R)
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
