#
#  Initialise.jl
#  VertexModel
#
#  Created by Christopher Revell on 31/01/2021.
#
#
# Function to initialise initial system

module Initialise

# Julia packages
using SparseArrays
using StaticArrays
using DelimitedFiles
#using DrWatson

# Local modules
include("InitialHexagons.jl"); using .InitialHexagons
include("VertexModelContainers.jl"); using .VertexModelContainers


function initialise(initialSystem,realTimetMax,γ,λ,preferredArea,pressureExternal,dt,tStar,outputTotal,t1Threshold)

    # Derived parameters
    tMax               = realTimetMax/tStar # Non dimensionalised maximum system run time
    outputInterval     = tMax/outputTotal   # Time interval for storing system data (non dimensionalised)
    preferredPerimeter = -λ/(2*γ)           # Cell preferred perimeter

    # Initialise system matrices from function or file
    # A is an incidence matrix mapping edges to vertices. Rows => edges; columns => vertices.
    # B is an incidence matrix mapping cells to edges. Rows => cells; columns => edges. (Values +/-1 for orientation)
    # C is an adjacency matrix (derived from A and B) mapping cells to vertices. Rows => cells; Columns => vertices. C=0.5*B̄*Ā
    # R is a vector of vertex positions in 2D. Each position stored as a 2-component static vector.
    if initialSystem=="single"
        A,B,R = initialHexagons(1)
    elseif initialSystem=="triple"
        A,B,R = initialHexagons(3)
    else
        # Import system matrices from file
        A = sparse(readdlm("data/input/$(initialSystem)_A.txt",',',Int64,'\n'))
        B = sparse(readdlm("data/input/$(initialSystem)_B.txt",',',Int64,'\n'))
        R0 = readdlm("data/input/$(initialSystem)_R.txt",',',Float64,'\n')
        R = Array{SVector{2,Float64}}(undef,size(A)[2])
        for i=1:size(R0)[1]
            R[i] = SVector{2}(R0[i,:])
        end
    end

    # Infer system information from matrices
    nCells            = size(B)[1]
    nEdges            = size(A)[1]
    nVerts            = size(A)[2]

    # Preallocate empty arrays for additional system matrices. See VertexModelContainers for explanations.
    Aᵀ                = spzeros(Int64,nVerts,nEdges)
    Ā                 = spzeros(Int64,nEdges,nVerts)
    Āᵀ                = spzeros(Int64,nVerts,nEdges)
    Bᵀ                = spzeros(Int64,nEdges,nCells)
    B̄                 = spzeros(Int64,nCells,nEdges)
    B̄ᵀ                = spzeros(Int64,nEdges,nCells)
    C                 = spzeros(Int64,nCells,nVerts)
    tempR             = Array{SVector{2,Float64}}(undef,nVerts)
    ΔR                = Array{SVector{2,Float64}}(undef,nVerts)
    cellEdgeCount     = zeros(Int64,nCells)
    boundaryVertices  = zeros(Int64,nVerts)
    cellPositions     = Array{SVector{2,Float64}}(undef,nCells)
    cellPerimeters    = zeros(nCells)
    cellOrientedAreas = Array{SMatrix{2,2,Float64}}(undef,nCells)
    cellAreas         = zeros(nCells)
    cellTensions      = zeros(nCells)
    cellPressures     = zeros(nCells)
    edgeLengths       = zeros(nEdges)
    edgeTangents      = Array{SVector{2,Float64}}(undef,nEdges)
    edgeMidpoints     = Array{SVector{2,Float64}}(undef,nEdges)
    vertexEdges       = Array{SVector{3,Float64}}(undef,nVerts)
    vertexCells       = Array{SVector{3,Float64}}(undef,nVerts)
    F                 = Array{SVector{2,Float64}}(undef,nVerts)
    rkCoefficients    = @SMatrix [
        0.0 0.5 0.5 0.5
        1.0 2.0 2.0 1.0
    ]
    ϵ                 = @SMatrix [
        0.0 1.0
        -1.0 0.0
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
        edgeLengths,
        edgeTangents,
        edgeMidpoints,
        vertexEdges,
        vertexCells,
        F,
        ϵ,
        rkCoefficients
    )

    # Pack parameters into a struct for convenience
    params = ParametersContainer(nVerts,
        nCells,
        nEdges,
        γ,
        λ,
        preferredPerimeter,
        preferredArea,
        pressureExternal,
        dt,
        outputInterval,
        tStar,
        realTimetMax,
        tMax,
        t1Threshold
    )

    return params,matrices

end

export initialise

end
