#
#  VertexModelContainers.jl
#  VertexModel
#
#  Created by Christopher Revell on 05/12/2021.
#
#
# Structs to hold system parameters and matrices, allowing us avoid passing dozens of arguments to each function

module VertexModelContainers

using SparseArrays
using StaticArrays

mutable struct ParametersContainer
    initialSystem      ::String   # System used for initialising simulations
    nVerts             ::Int64    # Number of cells
    nCells             ::Int64    # Number of edges
    nEdges             ::Int64    # Number of vertices
    γ                  ::Float64  # Parameter in energy relaxation
    λ                  ::Float64  # Parameter in energy relaxation
    L₀                 ::Float64  # Cell preferred perimeter length L₀ = -λ/(2*γ)
    A₀                 ::Float64  # Cell preferred area
    pressureExternal   ::Float64  # External pressure applied uniformly to system
    dt                 ::Float64  # Non dimensionalised time step
    outputTotal        ::Int64    # Total number of data outputs
    outputInterval     ::Float64  # Non dimensionalised data output interval
    viscousTimeScale   ::Float64  # Relaxation rate, approx from Sarah's data.
    realTimetMax       ::Float64  # Dimensionalised run time in seconds
    tMax               ::Float64  # Non dimensionalised run time
    realCycleTime      ::Float64  # Cell cycle time in seconds
    nonDimCycleTime    ::Float64  # Non dimensionalised cell cycle time
    t1Threshold        ::Float64  # Length of edge below which a T1 transition occurs
end

mutable struct MatricesContainer
    R                ::Vector{SVector{2, Float64}}    # Vector of 2D static vectors for each vertex (x,y) position
    tempR            ::Vector{SVector{2, Float64}}    # Vector of 2D static vectors to hold intermediate positions of each vertex in Runge-Kutta integration
    ΔR               ::Vector{SVector{2, Float64}}    # Vector of 2D static vectors to hold overall change of position of each vertex in Runge-Kutta integration
    A                ::SparseMatrixCSC{Int64, Int64}  # Incidence matrix mapping edges to vertices. Rows => edges; columns => vertices.
    B                ::SparseMatrixCSC{Int64, Int64}  # Incidence matrix mapping cells to edges. Rows => cells; columns => edges. (Values +/-1 for orientation)
    Aᵀ               ::SparseMatrixCSC{Int64, Int64}  # Transpose of incidence matrix A
    Ā                ::SparseMatrixCSC{Int64, Int64}  # Undirected adjacency matrix mapping edges to vertices (abs.(A))
    Āᵀ               ::SparseMatrixCSC{Int64, Int64}  # Transpose of adjacency matric Ā
    Bᵀ               ::SparseMatrixCSC{Int64, Int64}  # Transpose of incidence matrix B
    B̄                ::SparseMatrixCSC{Int64, Int64}  # Undirected adjacency matrix mapping cells to edges (abs.(B))
    B̄ᵀ               ::SparseMatrixCSC{Int64, Int64}  # Transpose of adjacency matric B̄
    C                ::SparseMatrixCSC{Int64, Int64}  # Adjacency matrix mapping cells to vertices. Rows => cells; Columns => vertices. (C=0.5*B̄*Ā)
    cellEdgeCount    ::Vector{Int64}                  # Vector of number of edges around each cell
    boundaryVertices ::Vector{Int64}                  # Vector of 1s and 0s denoting vertices that lie on the system boundary
    cellPositions    ::Vector{SVector{2, Float64}}    # Vector of 2D static vectors for each cell centre of mass
    cellPerimeters   ::Vector{Float64}                # Vector of scalar cell perimeter lengths
    cellOrientedAreas::Vector{SMatrix{2, 2, Float64}} # Vector of oriented cell areas. Each row is a 2x2 antisymmetric static matrix of the form [0 A / -A 0] where A is the scalar cell area
    cellAreas        ::Vector{Float64}                # Vector of scalar cell areas
    cellTensions     ::Vector{Float64}                # Vector of boundary tensions for each cell
    cellPressures    ::Vector{Float64}                # Vector of internal pressures for each cell
    cellAges         ::Vector{Float64}                # Vector of cell ages
    edgeLengths      ::Vector{Float64}                # Vector of lengths for each edge in the system
    edgeTangents     ::Vector{SVector{2, Float64}}    # Vector of 2D static vectors containing edge length and direction as a 2D vector
    edgeMidpoints    ::Vector{SVector{2, Float64}}    # Vector of 2D static vectors containing edge midpoints as (x,y) positions
    F                ::Matrix{SVector{2, Float64}}    # Matrix of 2D static vectors containing force vectors acting on each vertex and cell
    externalF        ::Vector{SVector{2, Float64}}    # Vector of 2D static vectors containing total force applied to each vertex by external pressure
    totalF           ::Vector{SVector{2, Float64}}    # Vector of 2D static vectors containing resultant force vectors acting on each vertex
    ϵ                ::SMatrix{2, 2, Float64, 4}      # Antisymmetric rotation matrix
    rkCoefficients   ::SMatrix{2, 4, Float64, 8}      # Matrix containing coefficients of Runge-Kutta integration
end

export ParametersContainer,MatricesContainer

end
