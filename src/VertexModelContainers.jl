#
#  VertexModelContainers.jl
#  VertexModel
#
#  Created by Christopher Revell on 05/12/2021.
#
#
# Structs to hold system parameters and matrices

module VertexModelContainers

using SparseArrays
using StaticArrays

struct ParametersContainer
    nVerts             ::Int64
    nCells             ::Int64
    nEdges             ::Int64
    γ                  ::Float64
    λ                  ::Float64
    preferredPerimeter ::Float64
    preferredArea      ::Float64
    pressureExternal   ::Float64
    dt                 ::Float64
    outputInterval     ::Float64
    tStar              ::Float64
    realTimetMax       ::Float64
    tMax               ::Float64
    t1Threshold        ::Float64
end

struct MatricesContainer
    A                ::SparseMatrixCSC{Int64, Int64}
    B                ::SparseMatrixCSC{Int64, Int64}
    Aᵀ               ::SparseMatrixCSC{Int64, Int64}
    Ā                ::SparseMatrixCSC{Int64, Int64}
    Āᵀ               ::SparseMatrixCSC{Int64, Int64}
    Bᵀ               ::SparseMatrixCSC{Int64, Int64}
    B̄                ::SparseMatrixCSC{Int64, Int64}
    B̄ᵀ               ::SparseMatrixCSC{Int64, Int64}
    C                ::SparseMatrixCSC{Int64, Int64}
    cellEdgeCount    ::Vector{Int64}
    boundaryVertices ::Vector{Int64}
    cellPositions    ::Vector{SVector{2, Float64}}
    cellPerimeters   ::Vector{Float64}
    cellOrientedAreas::Vector{SMatrix{2, 2, Float64}}
    cellAreas        ::Vector{Float64}
    cellTensions     ::Vector{Float64}
    cellPressures    ::Vector{Float64}
    edgeLengths      ::Vector{Float64}
    edgeTangents     ::Vector{SVector{2, Float64}}
    edgeMidpoints    ::Vector{SVector{2, Float64}}
    vertexEdges      ::Vector{SVector{3, Float64}}
    vertexCells      ::Vector{SVector{3, Float64}}
    F                ::Vector{SVector{2, Float64}}
    ϵ                ::SMatrix{2, 2, Float64, 4}
end

export ParametersContainer,MatricesContainer

end
