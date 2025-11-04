#
#  Initialise.jl
#  VertexModel
#
#  Created by Christopher Revell on 31/01/2021.
#
#
# Function to initialise vertex model system matrices and derived parameters

module Initialise

# Julia packages
using SparseArrays
using StaticArrays
using JLD2
using UnPack
using FromFile
using DrWatson
using Random
using Distributions
using Dates
using CircularArrays

# Local modules
@from "InitialSystemLayout.jl" using InitialSystemLayout
@from "VertexModelContainers.jl" using VertexModelContainers
@from "TopologyChange.jl" using TopologyChange
@from "SpatialData.jl" using SpatialData

function initialise(initialSystem,realTimetMax,γ,L₀,A₀,pressureExternal,viscousTimeScale,outputTotal,t1Threshold,realCycleTime,peripheralTension,setRandomSeed; nRows=9,modelChoice="quadratic",
    vertexWeighting=0, stretchType="none", realStretchTime=0, λs=0, κ=1, maxCells=1)

    # Calculate derived parameters
    tMax = realTimetMax / viscousTimeScale  # Non dimensionalised maximum system run time
    tStretch= realStretchTime/viscousTimeScale
    outputInterval = tMax / (outputTotal)     # Time interval for storing system data (non dimensionalised)
    λ = -2.0 * L₀ * γ
    nonDimCycleTime = realCycleTime / viscousTimeScale # Non dimensionalised cell cycle time

    # Set random seed value and allocate random number generator
    # Random seed set from current unix time, 
    # unless non zero value of setRandomSeed is passed, in which case random seed is passed value of setRandomSeed
    seed = (setRandomSeed == 0 ? floor(Int64, datetime2unix(now())) : setRandomSeed)
    Random.seed!(seed)

    distLogNormal = LogNormal(0.0, 0.2) #Distribution for division

    # Initialise system matrices from function or file
    if initialSystem in ["one", "three", "seven", "three_uneq", "three_neq2", "seven_eq", "sixteen"]
        # Create matrices for one, three, or seven cells geometrically
        A, B, R = initialCellConfig(initialSystem)
        cellTimeToDivide = rand(Uniform(0.0, nonDimCycleTime), size(B, 1))  # Random initial cell ages
    elseif initialSystem == "hex"
        A, B, R = initialSystemLayout(nRows, x=0)
        cellTimeToDivide = rand(Uniform(0.0, nonDimCycleTime), size(B, 1))  # Random initial cell ages
    elseif initialSystem == "hex_smooth"
        A, B, R = initialSystemLayout(nRows, edgeCells=true)
        cellTimeToDivide = rand(Uniform(0.0, nonDimCycleTime), size(B, 1))  # Random initial cell ages
    else
        # Import system matrices from final state of previous run
        importedData = load("$initialSystem"; 
            typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer"=>ParametersContainer, 
            "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer"=>MatricesContainer))
        @unpack A,B = importedData["matrices"]
        cellTimeToDivide = rand(Uniform(0.0,nonDimCycleTime),size(B,1))
        R = importedData["R"]
    end



    nCells = size(B, 1)
    nEdges = size(A, 1)
    nVerts = size(A, 2)
    #cellTimeToDivide = rand(LogNormal(0.0, 0.2), size(B, 1)).*rand(LogNormal(0.0, 0.2), size(B, 1)).*rand(LogNormal(0.0, 0.2), size(B, 1)).*rand(LogNormal(0.0, 0.2), size(B, 1)).*nonDimCycleTime
    # Fill preallocated matrices into struct for convenience
    matrices = MatricesContainer(
        A                 = A,
        B                 = B,
        Aᵀ                = spzeros(Int64, nVerts, nEdges),
        Ā                 = spzeros(Int64, nEdges, nVerts),
        Āᵀ                = spzeros(Int64, nVerts, nEdges),
        Bᵀ                = spzeros(Int64, nEdges, nCells),
        B̄                 = spzeros(Int64, nCells, nEdges),
        B̄ᵀ                = spzeros(Int64, nEdges, nCells),
        C                 = spzeros(Int64, nCells, nVerts),
        cellEdgeCount     = zeros(Int64, nCells),
        cellVertexOrders  = fill(CircularVector(Int64[]), nCells),
        cellEdgeOrders    = fill(CircularVector(Int64[]), nCells),
        boundaryVertices  = zeros(Int64, nVerts),
        boundaryEdges     = zeros(Int64, nEdges),
        cellPositions     = fill(SVector{2,Float64}(zeros(2)), nCells),
        cellPerimeters    = zeros(nCells),
        cellOrientedAreas = fill(SMatrix{2,2,Float64}(zeros(2,2)), nCells),
        cellAreas         = zeros(nCells),
        cellA₀s           = fill(A₀, nCells),
        cellL₀s           = fill(L₀, nCells),
        cellTensions      = zeros(nCells),
        cellPressures     = zeros(nCells),
        cellTimeToDivide  = cellTimeToDivide,
        μ                 = ones(nCells),
        Γ                 = γ.*ones(nCells),
        edgeLengths       = zeros(nEdges),
        edgeTangents      = fill(SVector{2,Float64}(zeros(2)), nEdges),
        edgeMidpoints     = fill(SVector{2,Float64}(zeros(2)), nEdges),
        edgeMidpointLinks = spzeros(SVector{2,Float64}, nCells, nVerts),
        timeSinceT1       = ones(nEdges).*(nonDimCycleTime/100) ,
        vertexAreas       = ones(nVerts),
        F                 = spzeros(SVector{2,Float64}, nVerts, nCells),
        externalF         = fill(SVector{2,Float64}(zeros(2)), nVerts),
        totalF            = fill(SVector{2,Float64}(zeros(2)), nVerts),
        ϵ                 = SMatrix{2, 2, Float64}([
                                0.0 1.0
                                -1.0 0.0
                            ]),
        cellShapeTensor   = fill(SMatrix{2,2,Float64}(zeros(2,2)), nCells),
        R_membrane         = fill(SVector{2,Float64}(zeros(2)), nVerts),
        Rt                 = fill(SVector{2,Float64}(zeros(2)), nVerts),
        R_final            = fill(SVector{2,Float64}(zeros(2)), nVerts),
        cellLineage       = collect(1:nCells),
        cellGeneration    = ones(nCells),
        cellIndex         = collect(1:nCells),   
    )

    # Pack parameters into a struct for convenience
    params = ParametersContainer(
        initialSystem     = initialSystem,
        nCells            = nCells,
        nEdges            = nEdges,
        nVerts            = nVerts,
        γ                 = γ,
        λ                 = λ,
        L₀                = L₀,
        A₀                = A₀,
        pressureExternal  = pressureExternal,
        outputTotal       = outputTotal,
        outputInterval    = outputInterval,
        viscousTimeScale  = viscousTimeScale,
        realTimetMax      = realTimetMax,
        tMax              = tMax,
        realCycleTime     = realCycleTime,
        nonDimCycleTime   = nonDimCycleTime,
        t1Threshold       = t1Threshold,
        peripheralTension = peripheralTension,
        seed              = seed,
        distLogNormal     = distLogNormal,
        modelChoice       = modelChoice,
        vertexWeighting   = vertexWeighting,
        λs                = λs,
        stretchType       = stretchType,
        realStretchTime   = realStretchTime,
        tStretch          = tStretch,
        tMemChange        = 0,
        κ                 = κ,
        maxCells          = maxCells
    )

    matrices.R_membrane.=R
    # Initial evaluation of matrices based on system topology
    topologyChange!(matrices)
    spatialData!(R, params, matrices)

    # Convert vector of SVectors to flat vector of Float64
    u0 = Float64[]
    for r in R
        append!(u0, r)
    end

    return u0, params, matrices

end

export initialise

end
