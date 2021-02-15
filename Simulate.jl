#
#  Simulate.jl
#  VertexModelJL
#
#  Created by Christopher Revell on 31/01/2021.
#
#

module Simulate

# Julia packages
using LinearAlgebra
using DelimitedFiles
using Plots
#using Random
#using Distributions
#using Base.Threads
#using StaticArrays

# Local modules
using TopologyChange
using CreateRunDirectory
using GeometryTools
using CalculateForce
using SingleHexagon
using Parameters

@inline @views function simulate(initialSystem)

    # Parameters


    A = readdlm("input/$(initialSystem)_A.txt", ' ', Float64, '\n') # Incidence matrix. Rows => edges; columns => vertices.
    B = readdlm("input/$(initialSystem)_B.txt", ' ', Float64, '\n') # Incidence matrix. Rows => cells; columns => edges. Values +/-1 for orientation
    R = readdlm("input/$(initialSystem)_R.txt", ' ', Float64, '\n')[:,1:2] # Coordinates of vertices
    #A,B,R = singleHexagon()
    # Infer system information from matrices
    nCells = size(B)[1]  # Number of cells
    nEdges = size(A)[1]  # Number of edges
    nVerts = size(A)[2]  # Number of vertices
    R.+= rand(Float64,nVerts,2).*0.2.-0.1
    R.*=2.0

    Aᵀ = zeros(Float64,nVerts,nEdges)
    Ā  = zeros(Float64,nEdges,nVerts)
    Āᵀ = zeros(Float64,nVerts,nEdges)
    Bᵀ = zeros(Float64,size(B)[2],size(B)[1])
    B̄  = zeros(Float64,size(B))
    B̄ᵀ = zeros(Float64,size(Bᵀ))
    C  = zeros(Float64,size(B)[1],size(A)[2])        # C adjacency matrix. Rows => cells; Columns => vertices
    cellEdgeCount          = zeros(Int64,nCells,1)   # Number of edges around each cell found by summing columns of B̄
    boundaryVertices       = zeros(Int64,nVerts,1)
    cellPositions          = zeros(Float64,nCells,2) # Cell centre positions
    cellPerimeters         = zeros(Float64,nCells,1)
    cellAreas              = zeros(Float64,nCells,1)
    cellOrientedAreas      = zeros(Float64,nCells,2,2)
    cellTensions           = zeros(Float64,nCells,1)
    cellPressures          = zeros(Float64,nCells,1)
    cellEffectivePressures = zeros(Float64,nCells,1)
    edgeLengths            = zeros(Float64,nEdges,1)
    edgeTangents           = zeros(Float64,nEdges,2)
    edgeMidpoints          = zeros(Float64,nEdges,2)
    cellOnes               = ones(Float64,nCells)    # Useful array for reusing in calculations
    F                      = zeros(Float64,nVerts,2)

    folderName = createRunDirectory(nCells,nEdges,nVerts,A,B,R)

    t = 0.00000000001
    topologyChange!(A,Ā,Aᵀ,Āᵀ,B,B̄,Bᵀ,B̄ᵀ,C,cellEdgeCount,boundaryVertices,cellOnes)

    @gif for t=0:dt:tMax

        geometryTools!(A,Ā,B,B̄,C,R,nCells,nEdges,nVerts,cellPositions,cellEdgeCount,cellAreas,cellOrientedAreas,cellPerimeters,cellTensions,cellPressures,cellEffectivePressures,edgeLengths,edgeMidpoints,edgeTangents,gamma,preferredPerimeter)

        calculateForce!(F,A,Ā,B,B̄,cellPressures,cellTensions,edgeTangents,edgeLengths,nVerts,nCells,nEdges)

        R .+= F.*dt

        scatter(R[:,1],R[:,2],xlims=(-2,2),ylims=(-2,2),aspect_ratio=:equal)

    end


end

export simulate

end
