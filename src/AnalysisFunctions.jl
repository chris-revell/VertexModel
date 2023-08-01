#
#  AnalysisFunctions.jl
#  VertexModel
#
#  Created by Christopher Revell on 31/01/2021.
#
#

module AnalysisFunctions

# Julia packages
using DrWatson
using LinearAlgebra
using SparseArrays
using StaticArrays
using UnPack
using GeometryBasics
using Random
using FromFile
using Colors

# Local modules
@from "$(srcdir("OrderAroundCell.jl"))" using OrderAroundCell

function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end

function makeCellPolygons(R,params,matrices)
    cellPolygons = Vector{Point2f}[]
    for i=1:params.nCells
        orderedVertices, orderedEdges = orderAroundCell(matrices,i)
        push!(cellPolygons,Point2f.(R[orderedVertices]))
    end
    return cellPolygons
end

function makeCellLinks(params,matrices)
    @unpack B,cellPositions,edgeMidpoints = matrices
    @unpack nCells,nEdges = params
    onesVec = ones(1,nCells)
    boundaryEdges = abs.(onesVec*B)
    cᵖ = boundaryEdges'.*edgeMidpoints
    T = SVector{2,Float64}[]
    for j=1:nEdges
        Tⱼ = @SVector zeros(2)
        for i=1:nCells
            Tⱼ = Tⱼ + B[i,j]*(cellPositions[i].-cᵖ[j])
        end
        push!(T,Tⱼ)
    end
    return T
end

function makeLinkTriangles(R,params,matrices)
    @unpack A,B,C,boundaryVertices,cellPositions,edgeMidpoints = matrices
    @unpack nCells,nVerts = params
    onesVec = ones(1,nCells)
    linkTriangles = Vector{Point2f}[]
    boundaryEdges = abs.(onesVec*B)
    for k=1:nVerts
        if boundaryVertices[k] == 0
            # If this vertex is not at the system boundary, link triangle is easily formed from the positions of surrounding cells
            vertexCells = findall(x->x!=0,C[:,k])
            push!(linkTriangles, Point2f.(cellPositions[vertexCells]))
        else
            # If this vertex is at the system boundary, we must form a more complex kite from surrounding cell centres and midpoints of surrounding boundary edges
            vertexCells = findall(x->x!=0,C[:,k])
            vertexEdges = findall(x->x!=0,A[:,k])
            boundaryVertexEdges = intersect(vertexEdges,findall(x->x!=0,boundaryEdges[1,:]))
            kiteVertices = [edgeMidpoints[boundaryVertexEdges]; cellPositions[vertexCells]]
            push!(kiteVertices,R[k])
            com = sum(kiteVertices)./length(kiteVertices)
            angles = Float64[]
            for p=1:length(kiteVertices)
                angle = atan((kiteVertices[p].-com)...)
                push!(angles,angle)
            end
            kiteVertices .= kiteVertices[sortperm(angles)]
            push!(linkTriangles,Point2f.(kiteVertices))
        end
    end
    return linkTriangles
end

function makeEdgeTrapezia(R,params,matrices)
    @unpack A,B,cellPositions = matrices
    @unpack nEdges = params
    edgeTrapezia = Vector{Point2f}[]
    for j=1:nEdges
        edgeCells = findall(x->x!=0,B[:,j])
        edgeVertices = findall(x->x!=0,A[j,:])
        if length(edgeCells) > 1
            push!(edgeTrapezia,Point2f.([R[edgeVertices[1]],cellPositions[edgeCells[1]],R[edgeVertices[2]],cellPositions[edgeCells[2]]]))
        else
            push!(edgeTrapezia,Point2f.([R[edgeVertices[1]],cellPositions[edgeCells[1]],R[edgeVertices[2]]]))
        end
    end
    return edgeTrapezia
end

function makeEdgeMidpointPolygons(params,matrices)
    for i=1:params.nCells
        orderedVertices, orderedEdges = orderAroundCell(matrices,i)
        push!(edgeMidpointPolygons,Point2f.(matrices.edgeMidpoints[orderedEdges]))
    end
    return edgeMidpointPolygons
end

# {curlᶜb}ᵢ
# Calculate curl on each cell
function calculateCellCurls(R,params,matrices)
    @unpack C,cellPositions,edgeMidpoints,edgeTangents,cellAreas = matrices
    @unpack nCells = params
    cellCurls = Float64[]
    for c=1:nCells
        cellVertices = findall(x->x!=0,C[c,:])
        vertexAngles = zeros(size(cellVertices))
        for (k,v) in enumerate(cellVertices)
            vertexAngles[k] = atan((R[v].-cellPositions[c])...)
        end
        m = minimum(vertexAngles)
        vertexAngles .-= m
        cellVertices .= cellVertices[sortperm(vertexAngles)]
        cellEdges = findall(x->x!=0,B[c,:])
        edgeAngles = zeros(size(cellEdges))
        for (k,e) in enumerate(cellEdges)
            edgeAngles[k] = atan((edgeMidpoints[e].-cellPositions[c])...)
        end
        edgeAngles .+= (2π-m)
        edgeAngles .= edgeAngles.%(2π)
        cellEdges .= cellEdges[sortperm(edgeAngles)]
        h = @SVector [0.0,0.0]
        curlSum = 0
        for (i,e) in enumerate(cellEdges)
            h = h + ϵ*F[cellVertices[i],c]
            curlSum += B[c,e]*(h⋅edgeTangents[e])/cellAreas[c]
        end
        push!(cellCurls,curlSum)
    end
    return cellCurls
end

# {divᶜb}ᵢ
# Calculate div on each cell
function calculateCellDivs(R,params,matrices)
    @unpack B,C,F,cellPositions,edgeMidpoints,edgeTangents,cellAreas,ϵ = matrices
    @unpack nCells = params
    cellDivs = Float64[]
    for c=1:nCells
        cellVertices = findall(x->x!=0,C[c,:])
        vertexAngles = zeros(size(cellVertices))
        for (k,v) in enumerate(cellVertices)
            vertexAngles[k] = atan((R[v].-cellPositions[c])...)
        end
        m = minimum(vertexAngles)
        vertexAngles .-= m
        cellVertices .= cellVertices[sortperm(vertexAngles)]
        cellEdges = findall(x->x!=0,B[c,:])
        edgeAngles = zeros(size(cellEdges))
        for (k,e) in enumerate(cellEdges)
            edgeAngles[k] = atan((edgeMidpoints[e].-cellPositions[c])...)
        end
        edgeAngles .+= (2π-m)
        edgeAngles .= edgeAngles.%(2π)
        cellEdges .= cellEdges[sortperm(edgeAngles)]
        h = @SVector [0.0,0.0]
        divSum = 0
        for (i,e) in enumerate(cellEdges)
            h = h + ϵ*F[cellVertices[i],c]
            divSum -= B[c,e]*(h⋅(ϵ*edgeTangents[e]))/cellAreas[c]
        end
        # divSum *= (-0.5)
        push!(cellDivs,divSum)
    end
    return cellDivs
end

# {divᵛb}ₖ
# Calculate div at each vertex
function calculateVertexDivs(R,params,matrices,q,linkTriangleAreas)
    @unpack C,cellPositions,F,ϵ = matrices
    @unpack nVerts = params
    vertexDivs = Float64[]
    for k=1:nVerts
        divSum = 0
        vertexCells = findall(x->x!=0,C[:,k])
        cellAngles = zeros(length(vertexCells))
        for i=1:length(cellAngles)
            cellAngles[i] = atan((cellPositions[vertexCells[i]].-R[k])...)
        end
        vertexCells .= vertexCells[sortperm(cellAngles,rev=true)]
        for i in vertexCells
            divSum += ((ϵ*q[i,k])⋅(ϵ*F[k,i]))/linkTriangleAreas[k]
        end
        push!(vertexDivs,divSum)
    end
    return vertexDivs
end

# {CURLᵛb}ₖ
# Calculate curl at each vertex
function calculateVertexCurls(R,params,matrices,q,linkTriangleAreas)
    @unpack C,cellPositions,F,ϵ = matrices
    @unpack initialSystem,nVerts,nCells,nEdges,γ,λ,pressureExternal,outputTotal,outputInterval,viscousTimeScale,realTimetMax,tMax,realCycleTime,nonDimCycleTime,t1Threshold = params
    vertexCurls = Float64[]
    # Working around a given vertex, an h force space point from a cell is mapped to the next edge anticlockwise from the cell
    for k=1:nVerts
        curlSum = 0
        vertexCells = findall(x->x!=0,C[:,k])
        cellAngles = zeros(length(vertexCells))
        for i=1:length(cellAngles)
            cellAngles[i] = atan((cellPositions[vertexCells[i]].-R[k])...)
        end
        vertexCells .= vertexCells[sortperm(cellAngles,rev=true)]
        for i in vertexCells
            curlSum += (q[i,k]⋅(ϵ*F[k,i]))/linkTriangleAreas[k]
        end
        push!(vertexCurls,curlSum)
    end
    return vertexCurls
end

function makeCellVerticesDict(params,matrices)
    cellVerticesDict = Dict()
    for i=1:params.nCells
        cellVertices, cellEdges = orderAroundCell(matrices, i)
        # Store sorted cell vertices for this cell
        cellVerticesDict[i] = cellVertices
    end
    return cellVerticesDict
end

function edgeLinkMidpoints(R,params,matrices,trapeziumAreas,T)
    @unpack A,B,cellPositions,edgeTangents,edgeMidpoints,ϵ = matrices
    @unpack nEdges = params

    # Rotation matrix around vertices is the opposite of that around cells
    ϵₖ = -1*ϵ
    onesVec = ones(1,nCells)
    boundaryEdges = abs.(onesVec*B)
    cᵖ = boundaryEdges'.*edgeMidpoints

    intersections = SVector{2,Float64}[]
    for j=1:nEdges
        if boundaryEdges[j] == 0
            k = findall(x->x<0,A[j,:])[1]
            i = findall(x->x<0,B[:,j])[1]
            mⱼ = R[k] .+ ((cellPositions[i].-R[k])⋅(ϵₖ*T[j]))/(2.0*trapeziumAreas[j]).*edgeTangents[j]
            push!(intersections,mⱼ)
        else
            push!(intersections,cᵖ[j])
        end
    end
    return intersections
end

function calculateSpokes(R,params,matrices)
    @unpack C,cellPositions = matrices
    @unpack nVerts,nCells = params
    q = Matrix{SVector{2,Float64}}(undef,nCells,nVerts)
    for i=1:nCells
        for k=1:nVerts
            q[i,k] = abs(C[i,k]).*(R[k].-cellPositions[i])
        end
    end
    return q
end


function calculateVertexMidpointCurls(params,matrices,intersections,linkTriangleAreas,q)
    @unpack A,B,ϵ = matrices
    @unpack nVerts,nCells,nEdges = params

    vertexMidpointCurls = Float64[]
    # Working around a given vertex, an h force space point from a cell is mapped to the next edge anticlockwise from the cell
    for k=1:nVerts
        curlSum = 0
        for i=1:nCells
            for j=1:nEdges
                curlSum -= B[i,j]*A[j,k]*(q[i,k]⋅intersections[j])/linkTriangleAreas[k]
            end
        end
        push!(vertexMidpointCurls,curlSum)
    end
    return vertexMidpointCurls
end

function calculateVertexMidpointDivs(params,matrices,intersections,linkTriangleAreas,q)
    @unpack A,B,ϵ = matrices
    @unpack nVerts,nCells,nEdges = params

    vertexMidpointDivs = Float64[]
    # Working around a given vertex, an h force space point from a cell is mapped to the next edge anticlockwise from the cell
    for k=1:nVerts
        divSum = 0
        for i=1:nCells
            for j=1:nEdges
                divSum -= B[i,j]*A[j,k]*((ϵ*q[i,k])⋅intersections[j])/linkTriangleAreas[k]
            end
        end
        push!(vertexMidpointDivs,divSum)
    end
    return vertexMidpointDivs
end

function calculateCellMidpointDivs(params,matrices,intersections,q)
    @unpack B,cellAreas,edgeTangents = matrices
    @unpack nCells = params
    cellMidpointDivs = Float64[]
    for i=1:nCells
        orderedVertices, orderedEdges = orderAroundCell(matrices,i)       
        divSum = 0
        for j in orderedEdges
            divSum -= B[i,j]*(intersections[j]⋅(ϵ*edgeTangents[j]))/cellAreas[i]
        end        
        push!(cellMidpointDivs,divSum)
    end
    return cellMidpointDivs
end

function calculateCellMidpointCurls(params,matrices,intersections,q)
    @unpack B,cellAreas,edgeTangents = matrices
    @unpack nCells = params
    cellMidpointCurls = Float64[]
    for i=1:nCells
        orderedVertices, orderedEdges = orderAroundCell(matrices,i)
        curlSum = 0
        for j in orderedEdges
            curlSum += B[i,j]*(intersections[j]⋅edgeTangents[j])/cellAreas[i]
        end
        push!(cellMidpointCurls,curlSum)
    end
    return cellMidpointCurls
end

export getRandomColor
export makeCellPolygons
export makeCellLinks
export makeLinkTriangles
export makeEdgeTrapezia
export makeEdgeMidpointPolygons
export calculateCellCurls
export calculateCellDivs
export calculateVertexDivs
export calculateVertexCurls
export makeCellVerticesDict
export edgeLinkMidpoints
export calculateSpokes
export calculateVertexMidpointCurls
export calculateVertexMidpointDivs
export calculateCellMidpointDivs
export calculateCellMidpointCurls

end #end module 