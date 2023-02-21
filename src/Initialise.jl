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

# Local modules
@from "$(projectdir("src","InitialHexagons.jl"))" using InitialHexagons
@from "$(projectdir("src","VertexModelContainers.jl"))" using VertexModelContainers
@from "$(projectdir("src","TopologyChange.jl"))" using TopologyChange
@from "$(projectdir("src","SpatialData.jl"))" using SpatialData
@from "$(projectdir("src","CalculateForce.jl"))" using CalculateForce

function initialise(initialSystem,realTimetMax,γ,L₀,A₀,pressureExternal,dt,viscousTimeScale,outputTotal,t1Threshold,realCycleTime,peripheralTension)

    # Calculate derived parameters
    tMax               = realTimetMax/viscousTimeScale  # Non dimensionalised maximum system run time
    outputInterval     = tMax/outputTotal               # Time interval for storing system data (non dimensionalised)
    λ = -2.0*L₀*γ
    nonDimCycleTime    = realCycleTime/viscousTimeScale # Non dimensionalised cell cycle time

    # Initialise system matrices from function or file
    if initialSystem in ["one","three","seven"]
        # Create matrices for one, three, or seven cells geometrically
        A,B,R = initialHexagons(initialSystem)
    else
        # Import system matrices from final state of previous run
        importedArrays = load("$initialSystem/matricesFinal.jld2")
        @unpack A,B,R,cellAges = importedArrays["matrices"]
    end

    # Infer system information from matrices
    nCells            = size(B)[1]
    nEdges            = size(A)[1]
    nVerts            = size(A)[2]

    # Preallocate empty arrays for additional system matrices.
    # See VertexModelContainers for explanations.
    Aᵀ                = spzeros(Int64,nVerts,nEdges)
    Ā                 = spzeros(Int64,nEdges,nVerts)
    Āᵀ                = spzeros(Int64,nVerts,nEdges)
    Bᵀ                = spzeros(Int64,nEdges,nCells)
    B̄                 = spzeros(Int64,nCells,nEdges)
    B̄ᵀ                = spzeros(Int64,nEdges,nCells)
    C                 = spzeros(Int64,nCells,nVerts)
    tempR             = Array{SVector{2,Float64}}(undef,nVerts)
    fill!(tempR,@SVector zeros(2))
    ΔR                = Array{SVector{2,Float64}}(undef,nVerts)
    fill!(ΔR,@SVector zeros(2))
    cellEdgeCount     = zeros(Int64,nCells)
    boundaryVertices  = zeros(Int64,nVerts)
    boundaryEdges     = zeros(Int64,nEdges)
    cellPositions     = Array{SVector{2,Float64}}(undef,nCells)
    fill!(cellPositions,@SVector zeros(2))
    cellPerimeters    = zeros(nCells)
    cellOrientedAreas = Array{SMatrix{2,2,Float64}}(undef,nCells)
    fill!(cellOrientedAreas,@SMatrix zeros(2,2))
    cellAreas         = zeros(nCells)
    cellTensions      = zeros(nCells)
    cellPressures     = zeros(nCells)
    initialSystem in ["one","three","seven"] ? cellAges = rand(nCells).*nonDimCycleTime : nothing  # Random initial cell ages
    edgeLengths       = zeros(nEdges)
    edgeTangents      = Vector{SVector{2,Float64}}(undef,nEdges)
    fill!(edgeTangents,@SVector zeros(2))
    edgeMidpoints     = Vector{SVector{2,Float64}}(undef,nEdges)
    fill!(edgeMidpoints,@SVector zeros(2))
    F                 = Matrix{SVector{2,Float64}}(undef,nVerts,nCells)
    fill!(F,@SVector zeros(2))
    externalF         = Array{SVector{2,Float64}}(undef,nVerts)
    fill!(externalF,@SVector zeros(2))
    totalF         = Array{SVector{2,Float64}}(undef,nVerts)
    fill!(totalF,@SVector zeros(2))
    rkCoefficients    = @SMatrix [  # Coefficients for Runge-Kutta integration
        0.0 0.5 0.5 0.5
        1.0 2.0 2.0 1.0
    ]
    ϵ                 = @SMatrix [  # Clockwise rotation matrix setting orientation of cell faces
        0.0 1.0
        -1.0 0.0
    ]

    # Pack matrices into a struct for convenience
    matrices = MatricesContainer(
        R,
        tempR,
        ΔR,
        A,
        B,
        Aᵀ,
        Ā,
        Āᵀ,
        Bᵀ,
        B̄,
        B̄ᵀ,
        C,
        cellEdgeCount,
        boundaryVertices,
        boundaryEdges,
        cellPositions,
        cellPerimeters,
        cellOrientedAreas,
        cellAreas,
        cellTensions,
        cellPressures,
        cellAges,
        edgeLengths,
        edgeTangents,
        edgeMidpoints,
        F,
        externalF,
        totalF,
        ϵ,
        rkCoefficients
    )

    # Pack parameters into a struct for convenience
    params = ParametersContainer(initialSystem,
        nVerts,
        nCells,
        nEdges,
        γ,
        λ,
        L₀,
        A₀,
        pressureExternal,
        dt,
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
    calculateForce!(R,params,matrices)

    return params,matrices

end

export initialise

end
