#
#  Initialise.jl
#  VertexModel
#
#  Function to initialise vertex model system matrices and derived parameters

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
# @from "InitialHexagons.jl" using InitialHexagons
@from "InitialSystemLayout.jl" using InitialSystemLayout
@from "VertexModelContainers.jl" using VertexModelContainers
@from "TopologyChange.jl" using TopologyChange
@from "SpatialData.jl" using SpatialData

function initialise(; initialSystem = "new",
        nCycles = 1,
        realCycleTime = 86400.0,
        realTimetMax = nCycles*realCycleTime,
        γ = 0.2,
        L₀ = 0.75,
        A₀ = 1.0,
        pressureExternal = 0.0,
        viscousTimeScale = 1000.0,
        outputTotal = 100,
        t1Threshold = 0.05,
        peripheralTension = 0.0,
        randomSeed = 0,
        nRows = 9,
        energyModel = "log",
        vertexWeighting = 1,
        R_in= spzeros(2),
        A_in= spzeros(2),
        B_in= spzeros(2),
        surfaceRadius = 20.0,
        surfaceReturnAmplitude = 100.0,
        initialEdgeLength = 0.75,
    )

    # Calculate derived parameters
    tMax = realTimetMax / viscousTimeScale  # Non dimensionalised maximum system run time
    outputInterval = tMax / (outputTotal-1)     # Time interval for storing system data (non dimensionalised)
    λ = -2.0 * L₀ * γ
    nonDimCycleTime = realCycleTime / viscousTimeScale # Non dimensionalised cell cycle time

    # Set random seed value and allocate random number generator
    # Random seed set from current unix time, 
    # unless non zero value of randomSeed is passed, in which case random seed is passed value of randomSeed
    seed = (randomSeed == 0 ? floor(Int64, datetime2unix(now())) : randomSeed)
    Random.seed!(seed)
    rng = MersenneTwister(seed)

    # Initialise system matrices from function or file
    if initialSystem == "new"
        isodd(nRows) && (nRows>1)  ? nothing : throw("nRows must be an odd number greater than 1.")
        A, B, R = initialSystemLayout(nRows=nRows, initialEdgeLength=initialEdgeLength=5*L₀/6)
        cellTimeToDivide = rand(rng,Uniform(0.0, nonDimCycleTime), size(B, 1))  # Random initial cell ages
    elseif initialSystem == "argument"
        R = R_in
        A = A_in
        B = B_in
        cellTimeToDivide = rand(rng,Uniform(0.0, nonDimCycleTime), size(B, 1))  # Random initial cell ages
    else
        # Import system matrices from final state of previous run
        importedData = load("$initialSystem"; 
            typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer"=>ParametersContainer, 
            "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer"=>MatricesContainer))
        @unpack A,B = importedData["matrices"]
        cellTimeToDivide = rand(rng,Uniform(0.0,nonDimCycleTime),size(B,1))
        R = importedData["R"]
    end

    nCells = size(B, 1)
    nEdges = size(A, 1)
    nVerts = size(A, 2)

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
        cellPositions     = fill(SVector{3,Float64}(zeros(3)), nCells),
        cellPerimeters    = zeros(nCells),
        # cellOrientedAreas = fill(SMatrix{2,2,Float64}(zeros(2,2)), nCells),
        cellAreas         = zeros(nCells),
        cellA₀s           = fill(A₀, nCells),
        cellL₀s           = fill(L₀, nCells),
        cellTensions      = zeros(nCells),
        cellPressures     = zeros(nCells),
        cellTimeToDivide  = cellTimeToDivide,
        μ                 = ones(nCells),
        Γ                 = γ.*ones(nCells),
        edgeLengths       = zeros(nEdges),
        edgeTangents      = fill(SVector{3,Float64}(zeros(3)), nEdges),
        edgeMidpoints     = fill(SVector{3,Float64}(zeros(3)), nEdges),
        edgeϵs            = fill(MMatrix{3,3,Float64}(zeros(3,3)), nEdges),
        edgeMidpointLinks = spzeros(SVector{3,Float64}, nCells, nVerts),
        timeSinceT1       = zeros(nEdges),
        vertexAreas       = ones(nVerts),
        F                 = spzeros(SVector{3,Float64}, nVerts, nCells),
        externalF         = fill(SVector{3,Float64}(zeros(3)), nVerts),
        totalF            = fill(SVector{3,Float64}(zeros(3)), nVerts),
        # ϵ                 = SMatrix{2, 2, Float64}([
        #                         0.0 1.0
        #                         -1.0 0.0
        #                     ]),
        # cellShapeTensor   = fill(SMatrix{2,2,Float64}(zeros(2,2)), nCells),
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
        nCycles           = nCycles,
        nonDimCycleTime   = nonDimCycleTime,
        t1Threshold       = t1Threshold,
        peripheralTension = peripheralTension,
        seed              = seed,
        rng               = rng,
        distLogNormal     = LogNormal(0.0, 0.2),
        energyModel       = energyModel,
        vertexWeighting   = vertexWeighting,
        surfaceCentre     = SVector{3,Float64}(0.0,0.0,-surfaceRadius),
        surfaceRadius     = surfaceRadius,
        surfaceReturnAmplitude = surfaceReturnAmplitude,
    )

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
