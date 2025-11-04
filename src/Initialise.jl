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
using CSV # For reading csv files
using DataFrames # For reading csv files
using MAT # For reading mat files

# Local modules
@from "initialSystemLayout.jl" using InitialSystemLayout
@from "VertexModelContainers.jl" using VertexModelContainers
@from "TopologyChange.jl" using TopologyChange
@from "SpatialData.jl" using SpatialData

# Extracting data from csv files: 
dfR = DataFrame(CSV.File(datadir("R.csv");header=false)) 
dfA = DataFrame(CSV.File(datadir("A.csv");header=false)) 
dfB = DataFrame(CSV.File(datadir("B.csv");header=false)) 
dfC = DataFrame(CSV.File(datadir("C.csv");header=false)) 


# Import the boundary cells: 
dfboundary = DataFrame(CSV.File(datadir("boundary_cells_indices.csv");header=false)) 
# Import the area vector: 
dfAreaVec = DataFrame(CSV.File(datadir("Area_vec.csv");header=false))
# Import cell centres: 
dfCentres = DataFrame(CSV.File(datadir("cellPositions.csv");header=false))
# Import edge vectors: 
dfEdges = DataFrame(CSV.File(datadir("edgeVectors.csv");header=false))
# Import link vectors: 
dfLinks = DataFrame(CSV.File(datadir("linkVectors.csv");header=false))

periodicR = [SVector(row[1],row[2]) for row in eachrow(dfR)]
dense_matrix_A = Matrix{Int}(dfA)
dense_matrix_B = Matrix{Int}(dfB)
dense_matrix_C = Matrix{Int}(dfC)
periodicA = sparse(dense_matrix_A)
periodicB = sparse(dense_matrix_B)
periodicBoundaryCellIndices = [Int(col[1]) for col in eachcol(dfboundary)]
cellPositionsPeriodic = [(Float64(row[1]), Float64(row[2])) for row in eachrow(dfCentres)]

function initialise(; initialSystem = "periodic",
        nCycles = 1,
        realCycleTime = 86400.0,
        realTimetMax = nCycles*realCycleTime,
        γ = 0.2,
        L0_A = 0.5,
        L0_B = 1.0,
        L₀ = 0.75,
        A₀ = 1.0,
        pressureExternal = 0.0,
        viscousTimeScale = 1000.0,
        outputTotal = 100,
        t1Threshold = 0.05,
        peripheralTension = 0.0,
        β,
        randomSeed = 0,
        nRows = 9,
        energyModel = "quadratic",
        vertexWeighting = 1,
        R_in= spzeros(2),
        A_in= spzeros(2),
        B_in= spzeros(2),
        L_x = 10,
        L_y = 10,
    )

    # Calculate derived parameters
    tMax = realTimetMax / viscousTimeScale  # Non dimensionalised maximum system run time
    outputInterval = tMax / (outputTotal-1)     # Time interval for storing system data (non dimensionalised)
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
        A, B, R = initialSystemLayout(nRows)
        cellTimeToDivide = rand(rng,Uniform(0.0, nonDimCycleTime), size(B, 1))  # Random initial cell ages
    elseif initialSystem == "periodic"
        A, B, R = periodicA, periodicB, periodicR
        cellTimeToDivide = rand(rng,Uniform(0.0, nonDimCycleTime), size(B, 1))  # Random initial cell ages
        cellPositions = cellPositionsPeriodic
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

    # Define the number of A and B cells
    nACells = Int(floor(nCells/2))

    cellsTypeA = randperm(rng, nCells)[1:nACells]   # random subset of cells
    cellsTypeB = setdiff(1:nCells, cellsTypeA)      # the remainder

    # Preferred perimeter depends on cell type
    # L0_A = 0.5    # preferred perimeter for Type A
    # L0_B = 1.0     # preferred perimeter for Type B

    cellL₀s = zeros(nCells)
    for i in 1:nCells
        if i in cellsTypeA
            cellL₀s[i] = L0_A
        else
            cellL₀s[i] = L0_B
        end
    end

    # Define contractility
    λs = -2.0 .* γ .* cellL₀s



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
        boundaryCells     = zeros(Int64, nCells),
        # cellPositions     = fill(SVector{2,Float64}(zeros(2)), nCells),
        cellPositions     = cellPositions,
        cellPerimeters    = zeros(nCells),
        cellOrientedAreas = fill(SMatrix{2,2,Float64}(zeros(2,2)), nCells),
        cellAreas         = zeros(nCells),
        cellA₀s           = fill(A₀, nCells),
        cellL₀s           = cellL₀s,
        cellTensions      = zeros(nCells),
        cellPressures     = zeros(nCells),
        cellTimeToDivide  = cellTimeToDivide,
        μ                 = ones(nCells),
        Γ                 = γ.*ones(nCells),
        edgeLengths       = zeros(nEdges),
        edgeTangents      = fill(SVector{2,Float64}(zeros(2)), nEdges),
        edgeMidpoints     = fill(SVector{2,Float64}(zeros(2)), nEdges),
        edgeMidpointLinks = spzeros(SVector{2,Float64}, nCells, nVerts),
        timeSinceT1       = zeros(nEdges),
        vertexAreas       = ones(nVerts),
        F                 = spzeros(SVector{2,Float64}, nVerts, nCells),
        externalF         = fill(SVector{2,Float64}(zeros(2)), nVerts),
        totalF            = fill(SVector{2,Float64}(zeros(2)), nVerts),
        ϵ                 = SMatrix{2, 2, Float64}([
                                0.0 1.0
                                -1.0 0.0
                            ]),
        cellShapeTensor   = fill(SMatrix{2,2,Float64}(zeros(2,2)), nCells),
    )

    # Pack parameters into a struct for convenience
    params = ParametersContainer(
        initialSystem     = initialSystem,
        nCells            = nCells,
        nEdges            = nEdges,
        nVerts            = nVerts,
        γ                 = γ,
        λs                = λs,
        L₀                = NaN,
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
        β                 = β,
        seed              = seed,
        rng               = rng,
        distLogNormal     = LogNormal(0.0, 0.2),
        energyModel       = energyModel,
        vertexWeighting   = vertexWeighting,
        cellsTypeA        = cellsTypeA,
        cellsTypeB        = cellsTypeB

    )

    # Initial evaluation of matrices based on system topology
    topologyChange!(R,params,matrices)
    spatialData!(R, params, matrices)

    # Checking the boundary cells match my Matlab vector: 
    println(length(periodicBoundaryCellIndices))
    println(count(!=(0), matrices.boundaryCells))
    # println(R)
    # println(matrices.cellPositions)
    

    # Convert vector of SVectors to flat vector of Float64
    u0 = Float64[]
    for r in R
        append!(u0, r)
    end

    return u0, params, matrices

end

export initialise

end
