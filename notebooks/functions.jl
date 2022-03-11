# Import Julia packages
using DrWatson
@quickactivate
using Revise
using LinearAlgebra
using SparseArrays
using StaticArrays
using UnPack
using GeometryBasics
using Random
using Colors
using JLD2
using Printf

# Local modules
includet("$(projectdir())/src/VertexModelContainers.jl"); using .VertexModelContainers


function cellPolygons(params,matrices)
    @unpack C,R,cellPositions = matrices
    @unpack nCells = params
    cellPolygons = Vector{Point2f}[]
    for i=1:nCells
        cellVertices = findall(x->x!=0,C[i,:])
        vertexAngles = zeros(size(cellVertices))
        for (k,v) in enumerate(cellVertices)
            vertexAngles[k] = atan((R[v].-cellPositions[i])...)
        end
        cellVertices .= cellVertices[sortperm(vertexAngles)]
        push!(cellPolygons,Point2f.(R[cellVertices]))
    end
    return cellPolygons
end

function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end
