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
@from "InitialHexagons.jl" using InitialHexagons
@from "largeInitialSystem.jl" using LargeInitialSystem
@from "VertexModelContainers.jl" using VertexModelContainers
@from "TopologyChange.jl" using TopologyChange
@from "SpatialData.jl" using SpatialData

function initialise(initialSystem,realTimetMax,γ,L₀,A₀,pressureExternal,viscousTimeScale,outputTotal,t1Threshold,realCycleTime,peripheralTension,setRandomSeed, surfaceRadius, surfaceReturnAmplitude)

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

    # Initialise system matrices from function or file
    if initialSystem in ["one", "three", "seven", "seven_original"]
        # Create matrices for one, three, or seven cells geometrically
        A, B, R = initialHexagons(initialSystem)
        cellTimeToDivide = rand(Uniform(0.0, nonDimCycleTime), size(B, 1))  # Random initial cell ages
    elseif initialSystem == "large"
        A, B, R = largeInitialSystem(2*L₀/6)
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

    # Fill preallocated matrices into struct for convenience
    matrices = MatricesContainer(
        A                   = A,
        B                   = B,
        Aᵀ                  = spzeros(Int64, nVerts, nEdges),                       # Aᵀ
        Ā                   = spzeros(Int64, nEdges, nVerts),                       # Ā
        Āᵀ                  = spzeros(Int64, nVerts, nEdges),                       # Āᵀ
        Bᵀ                  = spzeros(Int64, nEdges, nCells),                       # Bᵀ
        B̄                   = spzeros(Int64, nCells, nEdges),                       # B̄
        B̄ᵀ                  = spzeros(Int64, nEdges, nCells),                       # B̄ᵀ
        C                   = spzeros(Int64, nCells, nVerts),                       # C
        cellEdgeCount       = zeros(Int64, nCells),                                 # cellEdgeCount
        cellVertexOrders    = fill(CircularVector(Int64[]), nCells),                # cellVertexOrders
        cellEdgeOrders      = fill(CircularVector(Int64[]), nCells),                # cellEdgeOrders
        boundaryVertices    = zeros(Int64, nVerts),                                 # boundaryVertices
        boundaryEdges       = zeros(Int64, nEdges),                                 # boundaryEdges
        cellPositions       = fill(SVector{3,Float64}(zeros(3)), nCells),           # cellPositions
        cellPerimeters      = zeros(nCells),                                        # cellPerimeters
        # cellOrientedAreas   = fill(SMatrix{3,3,Float64}(zeros(3,3)), nCells),     # cellOrientedAreas
        cellϵs              = fill(SMatrix{3,3,Float64}(zeros(3,3)), nCells),       # cellOrientedAreas
        cellAreas           = zeros(nCells),                                        # cellAreas
        cellA₀s             = A₀.*ones(nCells),                                     # cellA₀s
        cellTensions        = zeros(nCells),                                        # cellTensions
        cellPressures       = zeros(nCells),                                        # cellPressures
        cellPerpAxes        = fill(SVector{3,Float64}(zeros(3)), nCells),           # cellPerpAxes
        cellTimeToDivide    = cellTimeToDivide,                                     # cellTimeToDivide
        μ                   = ones(nCells),                                         # μ
        Γ                   = γ.*ones(nCells),                                      # Γ
        edgeLengths         = zeros(nEdges),                                        # edgeLengths
        edgeTangents        = fill(SVector{3,Float64}(zeros(3)), nEdges),           # edgeTangents
        edgeMidpoints       = fill(SVector{3,Float64}(zeros(3)), nEdges),           # edgeMidpoints
        edgeMidpointLinks   = spzeros(SVector{3,Float64}, nCells, nVerts),          # edgeMidpointLinks
        timeSinceT1         = zeros(nEdges),                                        # timeSinceT1
        vertexAreas         = ones(nVerts),                                         # vertexAreas
        F                   = spzeros(SVector{3,Float64}, nVerts, nCells),          # F
        externalF           = fill(SVector{3,Float64}(zeros(3)), nVerts),           # externalF
        totalF              = fill(SVector{3,Float64}(zeros(3)), nVerts),           # totalF
        # cellShapeTensor     = fill(SMatrix{3,3,Float64}(zeros(3,3)), nCells),       # cellShapeTensor
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
        distLogNormal     = LogNormal(0.0, 0.2),
        surfaceCentre     = SVector{3,Float64}(0.0,0.0,-surfaceRadius),
        surfaceRadius     = surfaceRadius,
        surfaceReturnAmplitude = surfaceReturnAmplitude,
    )

    # Initial evaluation of matrices based on system topology
    topologyChange!(matrices)
    spatialData!(R, params, matrices)

    return R, params, matrices

end

export initialise

end
