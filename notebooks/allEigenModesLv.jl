using JLD2
using SparseArrays
using LinearAlgebra
using DrWatson
using FromFile
using UnPack
using CairoMakie
using Printf
using Colors
using CircularArrays
using GeometryBasics

@from "$(projectdir())/src/VertexModelContainers.jl" using VertexModelContainers
@from "$(projectdir())/src/AnalysisFunctions.jl" using AnalysisFunctions
@from "$(projectdir())/src/Laplacians.jl" using Laplacians
@from "$(projectdir())/src/Potentials.jl" using Potentials
@from "$(projectdir())/src/Eigenmodes.jl" using Eigenmodes
@from "$(projectdir())/src/OrderAroundCell.jl" using OrderAroundCell

frame = 100

folderName = "nCells=751_pressureExternal=0.5_realTimetMax=173000.0_stiffnessFactor=10.0_24-03-04-10-11-13"

files = [datadir("sims", folderName, "frameData", f) for f in readdir(datadir("sims", folderName, "frameData")) if occursin(".jld2", f)]

@unpack R, matrices, params = load(files[end];
    typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer" => ParametersContainer,
        "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer" => MatricesContainer))
@unpack B, Bᵀ, C, cellPositions = matrices
@unpack nCells, nVerts = params

fig = CairoMakie.Figure(size=(1000, 1000))
ax = Axis(fig[1, 1][1, 1], aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)

mkpath(datadir("sims", folderName,"eigenmodesLv","frame$(@sprintf("%03d", frame))"))

# cellPolygons = makeCellPolygonsOld(R,params,matrices)

# @unpack A,C,boundaryVertices,boundaryEdges,cellPositions,edgeMidpoints = matrices
# @unpack nVerts = params
# linkTriangles = Vector{Point{2,Float64}}[]
# cellEdgeOrders = fill(CircularVector(Int64[]),nCells)
# for i=1:length(cellEdgeOrders)
#     cellEdgeOrders[i] = orderAroundCell(matrices,i)[2]
# end
# for k=1:nVerts
#     if boundaryVertices[k] == 0
#         # If this vertex is not at the system boundary, link triangle is easily formed from the positions of surrounding cells
#         vertexCells = findall(x->x!=0,C[:,k])
#         push!(linkTriangles, Point{2,Float64}.(cellPositions[vertexCells]))
#     else
#         # If this vertex is at the system boundary, we must form a more complex kite from surrounding cell centres and midpoints of surrounding boundary edges
#         vertexCells = findall(x->x!=0, C[:,k])
#         vertexEdges = findall(x->x!=0, A[:,k])
#         boundaryVertexEdges = [v for v in vertexEdges if boundaryEdges[v]!=0] #vertexEdges ∩ findall(x->x!=0,boundaryEdges))
#         if length(vertexCells)>1
#             edge1 = (boundaryVertexEdges ∩ cellEdgeOrders[vertexCells[1]])[1]
#             edge2 = (boundaryVertexEdges ∩ cellEdgeOrders[vertexCells[2]])[1]
#             kiteVertices = [R[k], edgeMidpoints[edge1], cellPositions[vertexCells[1]], cellPositions[vertexCells[2]], edgeMidpoints[edge2]]
#         else 
#             kiteVertices = [R[k], edgeMidpoints[boundaryVertexEdges[1]], cellPositions[vertexCells[1]], edgeMidpoints[boundaryVertexEdges[2]]]
#         end
#         push!(linkTriangles,Point{2,Float64}.(kiteVertices))
#     end
# end

# edgeTrapezia = makeEdgeTrapezia(R,params,matrices)
# trapeziumAreas = abs.(area.(edgeTrapezia))
# linkTriangleAreas = abs.(area.(linkTriangles))    
# Lᵥ = makeLv(params,matrices,linkTriangleAreas,trapeziumAreas)
# decomposition = (eigen(Matrix(Lᵥ))).vectors
# mkpath(datadir("sims", folderName, "eigenmodesLv", "frame$(@sprintf("%03d", frame))"))

cellPolygons = makeCellPolygonsOld(R, params, matrices)
linkTriangles = makeLinkTriangles(R, params, matrices)

decomposition = eigenmodesLv(R, matrices, params)

for mode = 1:nVerts
    empty!(ax)
    lims = (-maximum(abs.(decomposition[:, mode])), maximum(abs.(decomposition[:, mode])))
    for k = 1:nVerts
        poly!(ax,
            linkTriangles[k],
            color=decomposition[k, mode],
            colorrange=lims,
            colormap=:bwr,
            strokewidth=1,
            strokecolor=(:black, 0.25),
        )
    end
    for i = 1:nCells
        poly!(ax,
            cellPolygons[i],
            color=(:white, 0.0),
            strokecolor=(:black, 1.0),
            strokewidth=1,
        )
    end
    save(datadir("sims", folderName, "eigenmodesLv", "frame$(@sprintf("%03d", frame))", "mode$(@sprintf("%03d", mode)).png"), fig)
end
