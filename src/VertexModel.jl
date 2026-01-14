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
@from "AnalysisFunctions.jl" using AnalysisFunctions




function vertexModel(;
    initialSystem = "periodic",
    cellLayout = "random",
    nRows = 9,
    nCycles = 0.01,
    realCycleTime = 86400.0,
    realTimetMax = nCycles*realCycleTime,
    γ = 0.04,
    L₀ = 0.5,
    A₀ = 1.0,
    viscousTimeScale = 1.0,
    pressureExternal = 0.0,
    peripheralTension = 0.0,
    t1Threshold = 0.05,
    β = 0.0,
    divisionToggle = 0,
    solver = SRIW1(),
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
    scatterVertices = 1,
    scatterCells = 1,
    plotForces = 0,
    plotEdgeMidpointLinks = 0,
    randomSeed = 0,
    abstol = 1e-6, 
    reltol = 1e-3,
    energyModel = "quadratic2pops",
    vertexWeighting = 1,
    R_in = spzeros(2),
    A_in = spzeros(2),
    B_in = spzeros(2), 
    L_x = 4,
    L_y = 4,
    Λ_00 = -0.04, #L0 = 0.5
    Λ_01 = 0.0,
    Λ_11 = -0.08, #L0 = 1.0
) # All arguments are optional and will be instantiated with these default values if not provided at runtime

    BLAS.set_num_threads(nBlasThreads)

    # Set up initial system, packaging parameters and matrices for system into params and matrices containers from VertexModelContainers.jl
    u0, params, matrices = initialise(initialSystem = initialSystem,
        cellLayout = cellLayout,
        nCycles = nCycles,
        realCycleTime = realCycleTime,
        realTimetMax = realTimetMax,
        γ = γ,
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

    ########################################################################################################################
    #           ADDING A GLOBAL TRY SO THAT MOVIE STILL GETS SAVED IF SOMETHING FAILS
    ########################################################################################################################
   
    try 

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

                energyComponent1 = sum((1/2).*(matrices.cellAreas .- 1).^2)
                energyComponent2 = sum((params.γ/2).*(matrices.cellPerimeters).^2)
                energyComponent3 = sum(matrices.Λs .* matrices.edgeLengths)

                # Print the three components of energy: 
                # println("Energy Components:", energyComponent1 , ",",energyComponent2, ",",energyComponent3, ", Total=", energyComponent1 .+ energyComponent2 .+ energyComponent3)
                
                
                
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
    
    catch err
        @error "VertexModel crashed — writing partial movie." exception=(err, catch_backtrace())

    finally

        # This will save the video even if an error occurs during the simulation

        if outputToggle == 1 && videoToggle == 1
            try
                save(datadir(folderName, "$(splitpath(folderName)[end]).mp4"), mov)
                @warn "Movie saved successfully (partial or complete)."
            catch saveErr
                @error "Movie failed to save in finally block." exception=(saveErr, catch_backtrace())
            end
        end

    end

    # If outputToggle==1, save animation object and save final system matrices
    (outputToggle == 1 && videoToggle == 1) ? save(datadir(folderName, "$(splitpath(folderName)[end]).mp4"), mov) : nothing

    P_eff = matrices.cellPressures .+ matrices.cellTensions.*matrices.cellPerimeters./(2.0.*matrices.cellAreas)
    println("Sum of effective cell pressures = ", sum(matrices.cellAreas.*P_eff))
    

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
