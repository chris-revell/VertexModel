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

# Local modules
@from "InitialHexagons.jl" using InitialHexagons
@from "largeInitialSystem.jl" using LargeInitialSystem
@from "VertexModelContainers.jl" using VertexModelContainers
@from "TopologyChange.jl" using TopologyChange
@from "SpatialData.jl" using SpatialData

function initialise(initialSystem,realTimetMax,γ,L₀,A₀,pressureExternal,viscousTimeScale,outputTotal,t1Threshold,realCycleTime,peripheralTension)

    # Calculate derived parameters
    tMax            = realTimetMax/viscousTimeScale  # Non dimensionalised maximum system run time
    outputInterval  = tMax/outputTotal               # Time interval for storing system data (non dimensionalised)
    λ               = -2.0*L₀*γ
    nonDimCycleTime = realCycleTime/viscousTimeScale # Non dimensionalised cell cycle time

    # Use this line if you want to force an identical random number sequences
    # rng = MersenneTwister(1234)

    # Initialise system matrices from function or file
    if initialSystem in ["one","three","seven"]
        # Create matrices for one, three, or seven cells geometrically
        A,B,R = initialHexagons(initialSystem)
        cellAges = rand(size(B,1)).*nonDimCycleTime  # Random initial cell ages
    elseif initialSystem=="large"
        A,B,R = largeInitialSystem()
        cellAges = rand(size(B,1)).*nonDimCycleTime  # Random initial cell ages
    else
        # Import system matrices from final state of previous run
        importedArrays = load("$initialSystem/dataFinal.jld2")
        @unpack A,B,cellAges = importedArrays["matrices"]
        R = importedArrays["R"]
    end

    nCells = size(B,1)
    nEdges = size(A,1)
    nVerts = size(A,2)

    # Pack preallocated matrices into a struct for convenience
    matrices = MatricesContainer(
        A,
        B,
        spzeros(Int64,nVerts,nEdges),                         # Aᵀ
        spzeros(Int64,nEdges,nVerts),                         # Ā
        spzeros(Int64,nVerts,nEdges),                         # Āᵀ
        spzeros(Int64,nEdges,nCells),                         # Bᵀ
        spzeros(Int64,nCells,nEdges),                         # B̄
        spzeros(Int64,nEdges,nCells),                         # B̄ᵀ
        spzeros(Int64,nCells,nVerts),                         # C
        zeros(Int64,nCells),                                  # cellEdgeCount
        zeros(Int64,nVerts),                                  # boundaryVertices
        zeros(Int64,nEdges),                                  # boundaryEdges
        fill(SVector{2,Float64}(zeros(2)), nCells),           # cellPositions
        zeros(nCells),                                        # cellPerimeters
        fill(SMatrix{2,2,Float64}(zeros(2,2)), nCells),       # cellOrientedAreas
        zeros(nCells),                                        # cellAreas
        zeros(nCells),                                        # cellTensions
        zeros(nCells),                                        # cellPressures
        cellAges,
        zeros(nEdges),                                        # edgeLengths
        fill(SVector{2,Float64}(zeros(2)), nEdges),           # edgeTangents
        fill(SVector{2,Float64}(zeros(2)), nEdges),           # edgeMidpoints
        fill(SVector{2, Float64}(zeros(2)), (nCells, nVerts)),# edgeMidpointLinks
        zeros(nEdges),                                        # timeSinceT1
        zeros(nVerts),                                        # vertexAreas
        fill(SVector{2,Float64}(zeros(2)), (nVerts, nCells)), # F
        fill(SVector{2,Float64}(zeros(2)), nVerts),           # externalF
        fill(SVector{2,Float64}(zeros(2)), nVerts),           # totalF
        SMatrix{2, 2, Float64}([                              # ϵ Clockwise rotation matrix setting orientation of cell faces
        0.0 1.0
        -1.0 0.0
        ])
    )

    # Pack parameters into a struct for convenience
    params = ParametersContainer(
        initialSystem,
        nCells,
        nEdges,
        nVerts,
        γ,
        λ,
        L₀,
        A₀,
        pressureExternal,
        outputTotal,
        outputInterval,
        viscousTimeScale,
        realTimetMax,
        tMax,
        realCycleTime,
        nonDimCycleTime,
        t1Threshold,
        peripheralTension
    )

    # Initial evaluation of matrices based on system topology
    topologyChange!(matrices)
    spatialData!(R,params,matrices)

    return R, params, matrices

end

export initialise

end
