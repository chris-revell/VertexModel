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
using StochasticDiffEq
using LinearAlgebra
using JLD2
using SparseArrays
using StaticArrays
using CairoMakie
using Printf
using DifferentialEquations

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
@from "Energy.jl" using Energy




function vertexModel(;
    initialSystem = "periodic",
    cellLayout = "random",
    nRows = 9,
    nCycles = 0.0001,
    realCycleTime = 86400.0,
    realTimetMax = nCycles*realCycleTime,
    γ = 0.04,
    L0_A = 0.5,
    L0_B = 0.5,
    L₀ = 0.5,
    A₀ = 1.0,
    viscousTimeScale = 1000.0,
    pressureExternal = 0.0,
    peripheralTension = 0.0,
    t1Threshold = 0.05,
    β = 0.0,
    divisionToggle = 0,
    solver = SRIW1(),
    nBlasThreads = 1,
    subFolder = "",
    outputTotal = 200,
    outputToggle = 1,
    frameDataToggle = 1,
    frameImageToggle = 1,
    printToggle = 1,
    videoToggle = 1,
    plotCells = 1,
    scatterEdges = 0,
    scatterVertices = 1,
    scatterCells = 1,
    plotForces = 0,
    plotEdgeMidpointLinks = 0,
    randomSeed = 0,
    abstol = 1e-7, 
    reltol = 1e-4,
    energyModel = "quadratic",
    vertexWeighting = 1,
    R_in = spzeros(2),
    A_in = spzeros(2),
    B_in = spzeros(2), 
    L_x = 4,
    L_y = 4,
    Λ_00 = 0.5,
    Λ_01 = 0.5,
    Λ_11 = 0.5,
) # All arguments are optional and will be instantiated with these default values if not provided at runtime

    BLAS.set_num_threads(nBlasThreads)

    # Set up initial system, packaging parameters and matrices for system into params and matrices containers from VertexModelContainers.jl
    u0, params, matrices = initialise(initialSystem = initialSystem,
        cellLayout = cellLayout,
        nCycles = nCycles,
        realCycleTime = realCycleTime,
        realTimetMax = realTimetMax,
        γ = γ,
        L0_A = L0_A,
        L0_B = L0_B,
        L₀ = L₀,
        A₀ = A₀,
        pressureExternal = pressureExternal,
        viscousTimeScale = viscousTimeScale,
        outputTotal = outputTotal,
        t1Threshold = t1Threshold,
        peripheralTension = peripheralTension,
        β = β,
        randomSeed = randomSeed,
        nRows = nRows,
        energyModel = energyModel,
        vertexWeighting = vertexWeighting,
        R_in = R_in,
        A_in = A_in,
        B_in = B_in,
        L_x=L_x,
        L_y=L_y,
        Λ_00 = Λ_00,
        Λ_01 = Λ_01,
        Λ_11 = Λ_11,
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
    prob = SDEProblem(model!, g!, u0, (0.0, Inf), (params, matrices))
    alltStops = collect(0.0:params.outputInterval:params.tMax) # Time points that the solver will be forced to land at during integration
    integrator = init(prob, solver; tstops=alltStops, abstol=abstol, reltol=reltol, save_on=false, save_start=false, save_end=true,verbose=true)
    outputCounter = [1]

   
    # Iterate until integrator time reaches max system time 
    while integrator.t <= params.tMax && (integrator.sol.retcode == ReturnCode.Default || integrator.sol.retcode == ReturnCode.Success)
        
        
        # Reinterpret state vector as a vector of SVectors 
        R = reinterpret(SVector{2,Float64}, integrator.u)
        if any(!isfinite, integrator.u)
            @show integrator.t
            @show integrator.u
            error("NaN or Inf detected in integrator.u")
        end

        # Note that reinterpreting accesses the same underlying data, so changes to R will update integrator.u and vice versa 

        # Output data to file 
        if integrator.t == alltStops[outputCounter[1]]
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

        # let E = energy(params, matrices)
        #     println("t=$(integrator.t): energy BEFORE step = ", E)
        # end

        # Step integrator forwards in time to update vertex positions 
        step!(integrator)

        # println("t=$(integrator.t): energy AFTER step = ", energy(params, matrices))

        if initialSystem == "periodic"
            # Wrap vertices into the periodic domain
            R = reinterpret(SVector{2,Float64}, integrator.u)
            for k in 1:length(R)
                x = R[k][1]
                y=R[k][2]
                # wrap 
                x = mod(x,L_x)
                y = mod(y,L_y)

                # Write back to state vector
                integrator.u[2k-1] = x 
                integrator.u[2k] = y
            end
        end

        # Update spatial data (edge lengths, cell areas, etc.) following iteration of the integrator
        spatialData!(R, params, matrices)

        # println("t=$(integrator.t): energy AFTER spatialData = ", energy(params, matrices))

        # Check system for T1 transitions 
        if t1Transitions!(integrator, params, matrices) > 0
            u_modified!(integrator, true)
            # senseCheck(matrices.A, matrices.B; marker="T1") # Check for nonzero values in B*A indicating error in incidence matrices           
            topologyChange!(R,params,matrices) # Update system matrices after T1 transition
            spatialData!(R, params, matrices) # Update spatial data after T1 transition  

            # println("t=$(integrator.t): energy AFTER T1 = ", energy(params, matrices))
        end
        if divisionToggle==1
            if division!(integrator, params, matrices) > 0
                u_modified!(integrator, true)
                # senseCheck(matrices.A, matrices.B; marker="division") # Check for nonzero values in B*A indicating error in incidence matrices          
                topologyChange!(R,params,matrices) # Update system matrices after division 
                spatialData!(R, params, matrices) # Update spatial data after division 
            end
        end
        # Update cell ages with (variable) timestep used in integration step
        matrices.cellTimeToDivide .-= integrator.dt
        matrices.timeSinceT1 .+= integrator.dt

        # println("Area sum=", sum(matrices.cellAreas))
    end

    # If outputToggle==1, save animation object and save final system matrices
    (outputToggle == 1 && videoToggle == 1) ? save(datadir(folderName, "$(splitpath(folderName)[end]).mp4"), mov) : nothing

    

    
    return integrator, matrices
end

# Function to load previously saved simulation data 
function loadData(relativePath; outputNumber=100)
    data = load(projectdir(relativePath, "frameData", "systemData$(@sprintf("%03d", outputNumber)).jld2"))
    return data["R"], data["matrices"], data["params"]
end

# Ensure code is precompiled
# @compile_workload begin
#     vertexModel(nCycles=0.01, outputToggle=0, frameDataToggle=0, frameImageToggle=0, printToggle=0, videoToggle=0)
# end

export vertexModel
export loadData 

end
