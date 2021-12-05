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
    tMax               = realTimetMax/tStar       # Non dimensionalised maximum system run time
    outputInterval     = tMax/outputTotal         # Time interval for storing system data (non dimensionalised)
    preferredPerimeter = -λ/(2*γ)        # Cell preferred perimeter
    ϵ                  = SMatrix{2,2}([0.0 1.0
                                       -1.0 0.0]) # Antisymmetric rotation matrix

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
    nCells            = size(B)[1]                               # Number of cells
    nEdges            = size(A)[1]                               # Number of edges
    nVerts            = size(A)[2]                               # Number of vertices

    # Preallocate empty arrays for additional system matrices
    Aᵀ                = spzeros(Int64,nVerts,nEdges)             # Transpose of incidence matrix A
    Ā                 = spzeros(Int64,nEdges,nVerts)             # Undirected adjacency matrix from absolute values of incidence matrix A
    Āᵀ                = spzeros(Int64,nVerts,nEdges)             # Undirected adjacency matrix from absolute values of transpose of incidence matrix Aᵀ
    Bᵀ                = spzeros(Int64,nEdges,nCells)             # Transpose of incidence matrix B
    B̄                 = spzeros(Int64,nCells,nEdges)             # Undirected adjacency matrix from absolute values of incidence matrix B
    B̄ᵀ                = spzeros(Int64,nEdges,nCells)             # Undirected adjacency matrix from absolute values of transpose of incidence matrix Bᵀ
    C                 = spzeros(Int64,nCells,nVerts)             # C adjacency matrix. Rows => cells; Columns => vertices. = 0.5*B̄*Ā
    tempR             = Array{SVector{2,Float64}}(undef,nVerts)  # Vector to store temporary positions in Runge-Kutta integration as static arrays
    ΔR                = Array{SVector{2,Float64}}(undef,nVerts)  # Vector to store change in R during Runge-Kutta integration as static arrays
    cellEdgeCount     = zeros(Int64,nCells)                      # 1D matrix containing number of edges around each cell, found by summing columns of B̄
    boundaryVertices  = zeros(Int64,nVerts)                      # 1D matrix containing labels of vertices at system boundary
    cellPositions     = Array{SVector{2,Float64}}(undef,nCells)  # Vector of 2D cell centre positions as static arrays
    cellPerimeters    = zeros(nCells)                            # 1D matrix of scalar cell perimeter lengths
    cellOrientedAreas = Array{SMatrix{2,2,Float64}}(undef,nCells)# Vector of oriented cell areas. Each row is a 2x2 antisymmetric static matrix of the form [0 A / -A 0] where A is the scalar cell area
    cellAreas         = zeros(nCells)                            # 1D matrix of scalar cell areas
    cellTensions      = zeros(nCells)                            # 1D matrix of scalar cell tensions
    cellPressures     = zeros(nCells)                            # 1D matrix of scalar cell pressures
    edgeLengths       = zeros(nEdges)                            # 1D matrix of scalar edge lengths
    edgeTangents      = Array{SVector{2,Float64}}(undef,nEdges)  # 2D matrix of tangent vectors for each edge (magnitude = edge length)
    edgeMidpoints     = Array{SVector{2,Float64}}(undef,nEdges)  # 2D matrix of position coordinates for each edge midpoint
    vertexEdges       = Array{SVector{3,Float64}}(undef,nVerts)  # 2D matrix containing the labels of all 3 edges around each vertex
    vertexCells       = Array{SVector{3,Float64}}(undef,nVerts)  # 2D matrix containing the labels of all 2-3 cells around each vertex
    F                 = Array{SVector{2,Float64}}(undef,nVerts)  # 3D array containing force vectors on vertex k from cell i, Fᵢₖ

    # Pack matrces into a struct for convenience
    matrices = MatricesContainer(A,B,Aᵀ,Ā,Āᵀ,Bᵀ,B̄,B̄ᵀ,C,cellEdgeCount,boundaryVertices,cellPositions,cellPerimeters,cellOrientedAreas,cellAreas,cellTensions,cellPressures,edgeLengths,edgeTangents,edgeMidpoints,vertexEdges,vertexCells,F,ϵ)
    
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

    return R,tempR,ΔR,params,matrices

end

export initialise

end
