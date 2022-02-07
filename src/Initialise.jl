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

# Local modules
include("InitialHexagons.jl"); using .InitialHexagons
include("VertexModelContainers.jl"); using .VertexModelContainers


function initialise(initialSystem,realTimetMax,γ,λ,preferredArea,pressureExternal,dt,viscousTimeScale,outputTotal,t1Threshold,realCycleTime)

    # Calculate derived parameters
    tMax               = realTimetMax/viscousTimeScale  # Non dimensionalised maximum system run time
    outputInterval     = tMax/outputTotal               # Time interval for storing system data (non dimensionalised)
    preferredPerimeter = -λ/(2*γ)                       # Preferred perimeter length for each cell
    nonDimCycleTime    = realCycleTime/viscousTimeScale # Non dimensionalised cell cycle time

    # Initialise system matrices from function or file
    if initialSystem in ["one","three","seven"]
        # Create matrices for one, three, or seven cells geometrically
        A,B,R = initialHexagons(initialSystem)
    else
        # Import system matrices from final state of previous run
        importedArrays = load("$initialSystem/matricesFinal.jld2")
        @unpack A,B,R = importedArrays["matrices"]
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
    cellPositions     = Array{SVector{2,Float64}}(undef,nCells)
    fill!(cellPositions,@SVector zeros(2))
    cellPerimeters    = zeros(nCells)
    cellOrientedAreas = Array{SMatrix{2,2,Float64}}(undef,nCells)
    fill!(cellOrientedAreas,@SMatrix zeros(2,2))
    cellAreas         = zeros(nCells)
    cellTensions      = zeros(nCells)
    cellPressures     = zeros(nCells)
    cellAges          = rand(nCells).*nonDimCycleTime   # Random initial cell ages
    cellEnergies      = zeros(nCells)
    edgeLengths       = zeros(nEdges)
    edgeTangents      = Vector{SVector{2,Float64}}(undef,nEdges)
    fill!(edgeTangents,@SVector zeros(2))
    edgeMidpoints     = Vector{SVector{2,Float64}}(undef,nEdges)
    fill!(edgeMidpoints,@SVector zeros(2))
    vertexEdges       = Array{SVector{3,Float64}}(undef,nVerts)
    fill!(vertexEdges,@SVector zeros(3))
    vertexCells       = Array{SVector{3,Float64}}(undef,nVerts)
    fill!(vertexCells,@SVector zeros(3))
    F                 = Matrix{SVector{2,Float64}}(undef,nVerts,nCells)
    fill!(F,@SVector zeros(2))
    externalF       = Array{SVector{2,Float64}}(undef,nVerts)
    fill!(externalF,@SVector zeros(2))
    rkCoefficients    = @SMatrix [  # Coefficients for Runge-Kutta integration
        0.0 0.5 0.5 0.5
        1.0 2.0 2.0 1.0
    ]
    ϵ                 = @SMatrix [  # Rotation matrix setting orientation of cell faces
        0.0 -1.0
        1.0 0.0
    ]

    # Pack matrces into a struct for convenience
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
        cellPositions,
        cellPerimeters,
        cellOrientedAreas,
        cellAreas,
        cellTensions,
        cellPressures,
        cellAges,
        cellEnergies,
        edgeLengths,
        edgeTangents,
        edgeMidpoints,
        vertexEdges,
        vertexCells,
        F,
        externalF,
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
        preferredPerimeter,
        preferredArea,
        pressureExternal,
        dt,
        outputTotal,
        outputInterval,
        viscousTimeScale,
        realTimetMax,
        tMax,
        realCycleTime,
        nonDimCycleTime,
        t1Threshold
    )

    return params,matrices

end

export initialise

end
