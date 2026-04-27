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
@from "InitialSystemLayoutPeriodic.jl" using InitialSystemLayoutPeriodic
@from "VertexModelContainers.jl" using VertexModelContainers
@from "TopologyChange.jl" using TopologyChange
@from "SpatialData.jl" using SpatialData

# println("mean area:",meanArea)

function initialise(; initialSystem,
        boundaryType,
        cellLayout,
        nCycles,
        realCycleTime,
        realTimetMax,
        γ,
        L₀,
        A₀ = 1.0,
        pressureExternal,
        viscousTimeScale,
        outputTotal,
        peripheralTension,
        β,
        randomSeed,
        nRows,
        energyModel,
        vertexWeighting,
        noiseWeighting,
        R_in,
        A_in,
        B_in,
        L_x,
        L_y,
        Λ_AA,
        Λ_AB,
        Λ_BB,
        Λ_AE,
        Λ_BE,
        Area_A_ratio,
        t1timeGap,
        spiky,
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
        A, B, R, t1Threshold, l_AA, l_BB, l_AB, l_AE, l_BE = initialSystemLayout(γ,Λ_AA,Λ_BB,Λ_AB,Λ_AE,Λ_BE, nRows,spiky)
        cellTimeToDivide = rand(rng,Uniform(0.0, nonDimCycleTime), size(B, 1))  # Random initial cell ages
        nCells = size(B, 1)
        nEdges = size(A, 1)
        nVerts = size(A, 2)

        println([l_AA, l_BB, l_AB, l_AE, l_BE])

        # Start by assigning all cells as type A:
        cellsTypeB = []
        cellsTypeA = []
        for i=1:nCells
            push!(cellsTypeA, i)
        end


        
    elseif initialSystem == "periodic"
        A,B,R,nACells,t1Threshold = initialSystemLayoutPeriodic(γ,L_x,L_y,Λ_AA,Λ_BB,Area_A_ratio)
        cellTimeToDivide = rand(rng,Uniform(0.0, nonDimCycleTime), size(B, 1))  # Random initial cell ages

        nCells = size(B, 1)
        println("N_c =" ,nCells)
        println("nACells = ", nACells)
        println("nBCells = ", nCells - nACells)
        nEdges = size(A, 1)
        nVerts = size(A, 2)


        cellsTypeA = 1:nACells 
        cellsTypeB = nACells+1:nCells 

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
        @unpack nCells,nEdges,nVerts,cellsTypeA,cellsTypeB,t1Threshold,Λ_AA,Λ_AB,Λ_BB,Λ_AE,Λ_BE,γ,l_AA,l_AB,l_BB,l_AE,l_BE = importedData["params"]
        

    end

    

    cellL₀s = L₀.*ones(nCells)
    

    # Label A cells as 0, B cells as 1.
    cellLabels = zeros(Int64, nCells)
    cellLabels[cellsTypeB] .= 1


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
        cellPositions     = fill(SVector{2,Float64}(zeros(2)), nCells),
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
                            ]), # clockwise orientation
        cellShapeTensor   = fill(SMatrix{2,2,Float64}(zeros(2,2)), nCells),
        cellLabels        = cellLabels,
        Λs                = zeros(nEdges),
        firstT1onEdge     = zeros(Int64, nEdges),
        P_effs            = zeros(nCells),
        T_effs            = zeros(nCells),
        ξs                = zeros(nCells),
        ξsVec             = fill(SVector{2,Float64}(zeros(2)), nCells),
        edgeLabels        = zeros(Int64, nEdges),
    )

    # Pack parameters into a struct for convenience
    params = ParametersContainer(
        initialSystem     = initialSystem,
        boundaryType      = boundaryType,
        cellLayout        = cellLayout,
        nCells            = nCells,
        nEdges            = nEdges,
        nVerts            = nVerts,
        γ                 = γ,
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
        β                 = β,
        seed              = seed,
        rng               = rng,
        distLogNormal     = LogNormal(0.0, 0.2),
        energyModel       = energyModel,
        vertexWeighting   = vertexWeighting,
        noiseWeighting    = noiseWeighting,
        cellsTypeA        = cellsTypeA,
        cellsTypeB        = cellsTypeB,
        L_x               = L_x,
        L_y               = L_y,
        Λ_AA              = Λ_AA,
        Λ_AB              = Λ_AB,
        Λ_BB              = Λ_BB,
        # Λ_AE              = 0.5*Λ_AA,
        # Λ_BE              = 0.5*Λ_BB,
        Λ_AE              = Λ_AE,
        Λ_BE              = Λ_BE,
        Area_A_ratio      = Area_A_ratio,
        t1timeGap         = t1timeGap,
        l_AA              = l_AA,
        l_BB              = l_BB,
        l_AB              = l_AB,
        l_AE              = l_AE,
        l_BE              = l_BE,
    )

    # Initial evaluation of matrices based on system topology
    topologyChange!(R,params,matrices)
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
