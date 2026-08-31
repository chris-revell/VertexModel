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
using Random
using DiscreteCalculus

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
@from "CubicSolutions.jl" using CubicSolutions
@from "ParameterDiagram.jl" using ParameterDiagram
@from "EdgeAblation.jl" using EdgeAblation





function vertexModel(;
    initialSystem = "new", # "new", "32-cell", "2-row" or jld2 path string 
    boundaryType = "free",
    cellLayout = "random",
    nRows = 3,
    nCycles = 1,
    realCycleTime = 400.0, # From Megan's data, using the division rate 0.15/min 
    realTimetMax = nCycles*realCycleTime,
    γ = 0.05,
    L₀ = 0.5,
    A₀ = 1.0,
    viscousTimeScale = 18.0,
    pressureExternal = 0.0,
    peripheralTension = 0.0,
    β = 0.0,
    divisionToggle = 0,
    ablationToggle = 0,
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
    scatterVertices = 0,
    scatterCells = 0,
    plotXis = 1,
    plotStresses = 1,
    plotForces = 0,
    plotEdgeMidpointLinks = 0,
    randomSeed = 0,
    energyModel = "quadratic2pops",
    vertexWeighting = 1,
    noiseWeighting = 1,
    R_in = spzeros(2),
    A_in = spzeros(2),
    B_in = spzeros(2), 
    L_x = 20,
    L_y = 20,
    Λ_AA = -0.2, 
    Λ_AB = -0.2,
    Λ_BB = -0.2,
    Λ_AE = -0.2, 
    Λ_BE = -0.2,
    Area_A_ratio = 0.5,
    t1timeGap = 5e-1,
    spiky = true,
    desiredNumCells = 100,
    plotOrientations = 1,
    edgeToAblate = [],
) # All arguments are optional and will be instantiated with these default values if not provided at runtime

    BLAS.set_num_threads(nBlasThreads)


    # Set up initial system, packaging parameters and matrices for system into params and matrices containers from VertexModelContainers.jl
    u0, params, matrices = initialise(initialSystem = initialSystem,
        boundaryType = boundaryType,
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
        peripheralTension = peripheralTension,
        β = β,
        randomSeed = randomSeed,
        nRows = nRows,
        energyModel = energyModel,
        vertexWeighting = vertexWeighting,
        noiseWeighting = noiseWeighting,
        R_in = R_in,
        A_in = A_in,
        B_in = B_in,
        L_x=L_x,
        L_y=L_y,
        Λ_AA = Λ_AA,
        Λ_AB = Λ_AB,
        Λ_BB = Λ_BB,
        Λ_AE = Λ_AE,
        Λ_BE = Λ_BE,
        Area_A_ratio = Area_A_ratio,
        t1timeGap = t1timeGap,
        spiky = spiky,
    )

    # Create directory in which to store date. Save parameters and store directory name for later use.
    if outputToggle == 1
        folderName = createRunDirectory(params,subFolder)
        # Create plot object for later use 
        if frameImageToggle==1 || videoToggle==1
            fig, ax1, ax2,ax3, cbar1,cbar2, mov = plotSetup()
        end
    end

    # Plot the parameter diagram: 
    fig2 = parameterDiagram(params)
    save(datadir(folderName, "parameterDiagram.png"), fig2)

    # Initialise plot of Peff, U, ξ:
    stressFig, axPeff, axU, axξ = stressPlotSetup()

    # Initialise a variable to store the energy at the previous step: 
    totalEnergyPrevious = 0.0

    # Flag to track whether the cell types have been assigned 
    cellsTypesAssigned = 0

    # Flag for tracking recoil after edge ablation
    trackInitialRecoil = 0

    # Initialise recoil plot: 
    recoilFig, axVelocity = recoilVelocityPlotSetup()

    
    # Flag to check if the first plot has been completed
    firstPlot = 0

    # Global try so that the movie still saves if there is an error:
    try 

        # Set up ODE integrator 
        # prob = SDEProblem(model!, g!, u0, (0.0, Inf), (params, matrices))
        # #alltStops = collect(0.0:params.outputInterval:params.tMax) # Time points that the solver will be forced to land at during integration
        # alltStops = collect(0.0:params.outputInterval:params.tMax)# Time points beyond which we plot the monolayer
        # integrator = init(prob, solver; abstol=abstol, reltol=reltol, save_on=false, save_start=false, save_end=true,verbose=true)
        
        if params.β ==0 
            abstol = 1e-9
            reltol = 1e-6
            solver = Tsit5()
            # Set ODE parameters: 
            prob = ODEProblem(model!,
                    u0,
                    (0.0, Inf),
                    (params, matrices)
                )
            alltStops = collect(0.0:params.outputInterval:params.tMax)# Time points beyond which we plot the monolayer
            integrator = init(prob,
                solver,
                tstops=alltStops,
                abstol=abstol,
                reltol=reltol,
                save_on=false,
                save_start=false,
                save_end=true,
            )  

            println("solver: Tsit5()")
        else
            abstol = 1e-6
            reltol = 1e-3
            solver = SRIW1()

            # Set up SDE integrator 
            prob = SDEProblem(model!, g!, u0, (0.0, Inf), (params, matrices))
            alltStops = collect(0.0:params.outputInterval:params.tMax)# Time points beyond which we plot the monolayer
            integrator = init(prob, solver; abstol=abstol, reltol=reltol, save_on=false, save_start=false, save_end=true,verbose=true)
        
            println("solver: SRIW1()")
        end
        
        outputCounter = [1]
        ablated = [false]
        
        # Iterate until integrator time reaches max system time 
        while integrator.t <= params.tMax && (integrator.sol.retcode == ReturnCode.Default || integrator.sol.retcode == ReturnCode.Success)
            if ablationToggle==1
                if outputCounter[1] > 33
                    break
                end
            end

            # Reinterpret state vector as a vector of SVectors 
            R = reinterpret(SVector{2,Float64}, integrator.u)
            if any(!isfinite, integrator.u)
                @show integrator.t
                @show integrator.u
                error("NaN or Inf detected in integrator.u")
            end

            # Note that reinterpreting accesses the same underlying data, so changes to R will update integrator.u and vice versa 
            
            # Output data to file 
            if integrator.t >= alltStops[outputCounter[1]] || (cellsTypesAssigned == 1 && initialSystem=="new") # To output the final state ones cell types are assigned
                
                if trackInitialRecoil ==1 

                    trackedDistance = trackVertices!(R,integrator.t,params,matrices)

                    recoilFig = recoilVelocityDiagram(integrator.t,trackedDistance,axVelocity,recoilFig)

                    save(datadir(folderName, "recoilVelocity.png"), recoilFig)
                end
                
                # Update progress on command line 
                printToggle == 1 ? println("$(@sprintf("%.2f", integrator.t))/$(@sprintf("%.2f", params.tMax)), $(outputCounter[1])/$outputTotal") : nothing            
                if frameDataToggle == 1
                    # Save system data to file 
                    jldsave(datadir(folderName, "frameData", "systemData$(@sprintf("%03d", outputCounter[1])).jld2"); matrices, params, R)
                end
                if frameImageToggle == 1 || videoToggle == 1
                    # Render visualisation of system and add frame to movie
                    cbar1, cbar2 = visualise(R, integrator.t, fig, ax1, ax2, ax3, cbar1, cbar2, mov, params, matrices, plotCells, scatterEdges, scatterVertices, scatterCells, plotXis, plotForces, plotEdgeMidpointLinks, plotStresses, trackInitialRecoil,plotOrientations)
                end
                # Save still image of this time step 
                frameImageToggle == 1 ? save(datadir(folderName, "frameImages", "frameImage$(@sprintf("%03d", outputCounter[1])).png"), fig) : nothing
                outputCounter[1] += 1

                # Plot the sum of P_effsA_i against time: 
                stressFig = P_effsDiagram(integrator.t,axPeff,stressFig,matrices)

                # Plot sum of U_i against time: 
                stressFig = U_iDiagram(integrator.t,axU,stressFig,matrices,params)

                # Plot sum of Aᵢξᵢ against time: 
                stressFig = ξ_iDiagram(integrator.t,axξ,stressFig,matrices,params)

                save(datadir(folderName, "stressPlots.png"), stressFig)

            end

            # Step integrator forwards in time to update vertex positions 
            step!(integrator)

            if boundaryType == "periodic"
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
                    topologyChange!(R,params,matrices) # Update system matrices after division 
                    spatialData!(R, params, matrices) # Update spatial data after division 
                end
            end
            # Update cell ages with (variable) timestep used in integration step
            matrices.cellTimeToDivide .-= integrator.dt
            matrices.timeSinceT1 .+= integrator.dt

            if ablationToggle == 1 && firstPlot == 1
                # if !(ablated[1]) && integrator.t>params.tMax/10.0
                if !(ablated[1])
                    systemCOM = sum(R)./params.nVerts 
                    # params.jAblated = findmin([norm(matrices.edgeMidpoints[j].-systemCOM) for j=1:params.nEdges])[2] # central edge
                    params.jAblated = edgeToAblate
                    println("jAblated=",params.jAblated)
                    # Find the vertices to track 
                    params.k_tracked = findall(x -> x!=0, @view matrices.A[params.jAblated,:])
                    println(params.k_tracked)
                    edgeAblation!(params.jAblated, params, matrices, integrator)
                    topologyChange!(R,params,matrices)
                    spatialData!(R, params, matrices) # Update spatial data after T1 transition  
                    ablated[1] = true

                    divisionToggle=0 # Stop divisions after ablation so we can see the effect clearly without the system getting too big
                    trackInitialRecoil = 1 # Track vertex positions for recoil velocity
                    
                end

                
            end

            
            
            

        

            # Check if we have assigned cell types in the previous run
            if cellsTypesAssigned ==1 && initialSystem == "new"
                break
            end 

            # For the case where we grow the monolayer: 
            if initialSystem == "new" && boundaryType == "free" && params.nCells >= desiredNumCells && cellsTypesAssigned == 0
                println("Desired cell number reached. Assigning cell types")

                @unpack nCells,Area_A_ratio = params
                nACells = floor(Int64, nCells*Area_A_ratio)
                params.cellsTypeB = []
                params.cellsTypeA = randperm(nCells)[1:nACells]
                for i=1:nCells
                    if i in params.cellsTypeA
                    else
                        push!(params.cellsTypeB, i)
                    end
                end
                matrices.cellLabels = zeros(Int64, nCells)
                matrices.cellLabels[params.cellsTypeB] .= 1
                cellsTypesAssigned = 1
            elseif initialSystem == "32-cell" && boundaryType == "free" && params.nCells >= 15 && cellsTypesAssigned ==0
                println("32 cells reached. Assign one type-B cell")
                @unpack nCells = params
                nACells = nCells - 1
                systemCOM = sum(matrices.cellPositions)./params.nCells 
                iSelected = findmin([norm(matrices.cellPositions[i].-systemCOM) for i=1:params.nCells])[2] # central edge
                params.cellsTypeB = [iSelected]
                params.cellsTypeA = setdiff(1:nCells, params.cellsTypeB)
                matrices.cellLabels = zeros(Int64, nCells)
                matrices.cellLabels[params.cellsTypeB] .= 1
                cellsTypesAssigned = 1
                # matrices.cellTimeToDivide[iSelected] = rand(params.distLogNormal)*params.nonDimCycleTime

            end

            firstPlot = 1

            
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

function ablationLoop(jld2PathString,edgeVector)

    for edge in edgeVector
        println("Ablating edge $edge")
        # vertexModel(initialSystem = jld2PathString,boundaryType="free",ablationToggle=1,β=0.0,nCycles=0.001,edgeToAblate=edge)
        vertexModel(initialSystem = jld2PathString,boundaryType="free",ablationToggle=1,β=0.0,nCycles=0.1,edgeToAblate=edge)
    end
    
end

function extractRecoilVecs(jld2PathStrings)

    DistVecs = Vector{Vector{Float64}}()
    timeVec = []
    # For error bars: 
    # maxDistVec = zeros(30)
    # minDistVec = zeros(30)
    firstVec = true 

    for jld2PathString in jld2PathStrings
        df = load(jld2PathString,"matrices")
        timeVecTemp = df.trackedTimePoints
        distVec = df.trackedVertDistance

        timeVecTemp = timeVecTemp[1:30]
        distVec = distVec[1:30]

        timeVecTemp = timeVecTemp .- timeVecTemp[1] # Starting from 0
        distVec = distVec .- distVec[1] # Starting from 0
        if firstVec
            timeVec = timeVecTemp
            firstVec = false
        end

        push!(DistVecs, distVec)
        println("distVec = ",distVec)
    end
    println("timeVec = ",timeVec)
    return nothing
end


function recoilComparisonPlot(timeVec,distVecs,lowerErrVecs,upperErrVecs,vecLabels)

    # Vec labels = 0,1,2 for edge type AA,BB,AB respectively. 

    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    recoilFig = Figure(size=(600,600))

    # Initialise a figure for tracking sum of P_effsA_i: 
    gridRecoil = recoilFig[1,1] = GridLayout()
    axVelocity = Axis(gridRecoil[1,1],aspect=1)
    axVelocity.title = "Initial recoil comparison averaged over 4 simulations"
    axVelocity.xlabel = "t"
    axVelocity.ylabel = "Δ|r₁-r₂|"

    for i in eachindex(distVecs)
        vec = distVecs[i]
        if vecLabels[i] == 0
            lines!(axVelocity,timeVec,vec,color=:blue, label="AA")
            errorbars!(axVelocity, timeVec, vec, lowerErrVecs[i], upperErrVecs[i], color=:blue)
        elseif vecLabels[i] == 1
            lines!(axVelocity,timeVec,vec,color=:orange, label="BB")
            errorbars!(axVelocity, timeVec, vec, lowerErrVecs[i], upperErrVecs[i], color=:orange)
        elseif vecLabels[i] == 2
            lines!(axVelocity,timeVec,vec,color=:grey, label="AB")
            # errorbars!(axVelocity, timeVec, vec, lowerErrVecs[i], upperErrVecs[i], color=:grey)
        end
    end
    axislegend(axVelocity)
    display(recoilFig)
    # save(datadir(folderName, "recoilComparison.png"), recoilFig)

end # end function 

function computeCoupleStressesFromSimulation(;jld2pathString, plotForces)
    # Function to work with DiscreteCalculus.jl couple stress functions, from a simulation in equilbrium. 
    # Insure jld2string input is NOT relative 
    R = load(jld2pathString,"R")
    params = load(jld2pathString,"params")
    matrices = load(jld2pathString,"matrices")
    @unpack A,B,C,F = matrices

    h = hNetwork(R,A,B,F)
    coupleStresses = -curlᵛ(R,A,B,h)
    vertexTriangles = findCellLinkVertexTriangles(R,A,B)

    # Plot couple stresses on vertices: 
    coupleStressFig, cellTypeAx, PeffAx, ξAx, coupleStressAx, PeffCbar, ξCbar, coupleStressCbar = coupleStressPlotSetup()

    PeffCbar, ξCbar, coupleStressCbar = visualiseCoupleStresses(R, coupleStressFig, cellTypeAx, coupleStressAx, coupleStressCbar, params, matrices, coupleStresses, vertexTriangles)

    # Find the folder the data has been taken from: 
    folderName = dirname(dirname(jld2pathString))
    mkpath(datadir(folderName))
    save(datadir(folderName, "coupleStressPlot.png"), coupleStressFig)

    if plotForces == 1
        # Case of plotting a cluster of cells to check the rotated forces align with the vector force potential; 
        forceComparisonFig, forceAx, polyAx, forceComparisonCbar = forceComparisonPlotSetup()
        # Select the cluster of cells we would like to consider: 
        iC = rand(1:params.nCells, 1)[1] 
        # iC = 2
        # Determine associated vertices: 
        kC = findall(x -> x!=0, @view C[iC,:])
        clusterCells = []
        # Determine this cell's neighbours: 
        for k in kC
            cells = findall(x -> x!=0, @view C[:,k])
            append!(clusterCells,cells)
        end
        unique!(clusterCells)

        forceComparisonCbar = visualiseForceComparison(R, h, clusterCells, forceComparisonFig, polyAx,forceAx, forceComparisonCbar, params, matrices, coupleStresses, vertexTriangles)
        save(datadir(folderName, "forceComparisonPlot.png"), forceComparisonFig)
    end

    

    return nothing
end

export vertexModel
export loadData 
export extractRecoilVecs, recoilComparisonPlot, ablationLoop, computeCoupleStressesFromSimulation

end
