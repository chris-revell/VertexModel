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

function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end

function makeCellPolygons(params,matrices)
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

function makeCellLinks(params,matrices)
    @unpack B, edgeMidpoints, cellPositions = matrices
    @unpack nCells = params
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

function makeLinkTriangles(params,matrices)
    @unpack boundaryVertices, C, cellPositions, A, B = matrices
    @unpack nVerts = params
    onesVec = ones(1,nCells)
    linkTriangles = Vector{Point2f}[]
    boundaryEdges = abs.(onesVec*B)
    for k=1:nVerts
        if boundaryVertices[k] == 0
            vertexCells = findall(x->x!=0,C[:,k])
            push!(linkTriangles, Point2f.(cellPositions[vertexCells]))
        else
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

function makeEdgeTrapezia(params,matrices)
    @unpack A, B, cellPositions = matrices
    @unpack nEdges = params
    edgeTrapezia = Vector{Point2f}[]
    for j=1:nEdges
        edgeCells = findall(x->x!=0,B[:,j])
        edgeVertices = findall(x->x!=0,A[j,:])
        trapeziumVertices = [R[edgeVertices]; cellPositions[edgeCells]]
        com = sum(trapeziumVertices)./length(trapeziumVertices)
        angles = Float64[]
        for p=1:length(trapeziumVertices)
            angle = atan((trapeziumVertices[p].-com)...)
            push!(angles,angle)
        end
        trapeziumVertices .= trapeziumVertices[sortperm(angles)]
        push!(edgeTrapezia,Point2f.(trapeziumVertices))
    end
    return edgeTrapezia
end

function makeEdgeMidpointPolygons(params,matrices)
    @unpack B, edgeMidpoints, cellPositions = matrices
    @unpack nCells = params
    edgeMidpointPolygons = Vector{Point2f}[]
    for i=1:nCells
        cellEdges = findall(!iszero,B[i,:])
        edgeAngles = zeros(size(cellEdges))
        for (k,v) in enumerate(cellEdges)
            edgeAngles[k] = atan((edgeMidpoints[v].-cellPositions[i])...)
        end
        cellEdges .= cellEdges[sortperm(edgeAngles)]
        push!(edgeMidpointPolygons,Point2f.(edgeMidpoints[cellEdges]))
    end
    return edgeMidpointPolygons
end


# {curlᶜb}ᵢ
# Calculate curl on each cell
function calculateCellCurls(params,matrices)
    @unpack B, C, R, F, cellPositions, edgeMidpoints, edgeTangents, cellAreas = matrices
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
function calculateCellDivs(params,matrices)
    @unpack B, C, R, ϵ, cellPositions, cellAreas = matrices
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
        divSum *= (-0.5)
        push!(cellDivs,divSum)
    end
    return cellDivs
end

# {divᵛb}ₖ
# Calculate div at each vertex
function calculateVertexDivs(params,matrices,T,linkTriangleAreas)
    @unpack A, C, R, edgeTangents, ϵ, boundaryVertices = matrices
    @unpack nVerts = params

    # Rotation matrix around vertices is the opposite of that around cells
    ϵₖ = -1*ϵ

    vertexDivs = Float64[]
    # Working around a given vertex, an h force space point from a cell is mapped to the next edge anticlockwise from the cell
    for k=1:nVerts
        if boundaryVertices[k] == 0
            vertexEdges = findall(x->x!=0,A[:,k])
            edgeAngles = zeros(length(vertexEdges))
            for (i,e) in enumerate(vertexEdges)
                edgeAngles[i] = atan((-A[e,k].*edgeTangents[e])...)
            end
            m = minimum(edgeAngles)
            edgeAngles .-= m
            vertexEdges .= vertexEdges[sortperm(edgeAngles,rev=true)]

            vertexCells = findall(x->x!=0,C[:,k])
            cellAngles = zeros(length(vertexCells))
            for i=1:length(cellAngles)
                cellAngles[i] = atan((cellPositions[vertexCells[i]].-R[k])...)
            end
            cellAngles .+= 2π-m
            cellAngles .= cellAngles.%2π
            vertexCells .= vertexCells[sortperm(cellAngles,rev=true)]
            h = @SVector [0.0,0.0]
            divSum = 0
            for (i,j) in enumerate(vertexEdges)
                h = h + ϵ*F[k,vertexCells[i]]
                divSum -= A[j,k]*((ϵₖ*T[j])⋅h)/linkTriangleAreas[k]
            end
        else
            # Set angles relative to pressure force angle, equivalent to the angle of a cell that doesn't actually exist
            # pressureAngle = atan((-1.0.*externalF[k])...)
            # vertexEdges = findall(x->x!=0,A[:,k])
            # edgeAngles = zeros(length(vertexEdges))
            # for (i,e) in enumerate(vertexEdges)
            #     edgeAngles[i] = atan((-A[e,k].*edgeTangents[e])...)
            # end
            #
            # edgeAngles .+= 2π-pressureAngle
            # edgeAngles .= edgeAngles.%2π
            # vertexEdges .= vertexEdges[sortperm(edgeAngles,rev=true)]
            #
            # vertexCells = findall(x->x!=0,C[:,k])
            # cellAngles = zeros(length(vertexCells))
            # for i=1:length(cellAngles)
            #     cellAngles[i] = atan((cellPositions[vertexCells[i]].-R[k])...)
            # end
            # cellAngles .+= 2π-pressureAngle
            # cellAngles .= cellAngles.%2π
            # vertexCells .= vertexCells[sortperm(cellAngles,rev=true)]
            #
            # h = @SVector [0.0,0.0]
            divSum = 0
            # h = h + ϵ*externalF[k]
            # divSum -= A[vertexEdges[1],k]*((ϵₖ*T[vertexEdges[1]])⋅h)/linkTriangleAreas[k]
            # for (i,j) in enumerate(vertexEdges[2:end])
            #     h = h + ϵ*F[k,vertexCells[i]]
            #     divSum -= A[j,k]*((ϵₖ*T[j])⋅h)/linkTriangleAreas[k]
            # end
        end
        divSum *= (-0.5)
        push!(vertexDivs,divSum)
    end
    return vertexDivs
end

# {CURLᵛb}ₖ
# Calculate curl at each vertex
function calculateVertexCurls(params,matrices,T,linkTriangleAreas)
    @unpack A, C, R, edgeTangents, ϵ, boundaryVertices = matrices
    @unpack nVerts = params

    vertexCurls = Float64[]
    # Working around a given vertex, an h force space point from a cell is mapped to the next edge anticlockwise from the cell
    for k=1:nVerts
        if boundaryVertices[k] == 0
            vertexEdges = findall(x->x!=0,A[:,k])
            edgeAngles = zeros(length(vertexEdges))
            for (i,e) in enumerate(vertexEdges)
                edgeAngles[i] = atan((-A[e,k].*edgeTangents[e])...)
            end
            m = minimum(edgeAngles)
            edgeAngles .-= m
            vertexEdges .= vertexEdges[sortperm(edgeAngles,rev=true)]
            vertexCells = findall(x->x!=0,C[:,k])
            cellAngles = zeros(length(vertexCells))
            for i=1:length(cellAngles)
                cellAngles[i] = atan((cellPositions[vertexCells[i]].-R[k])...)
            end
            cellAngles .+= 2π-m
            cellAngles .= cellAngles.%2π
            vertexCells .= vertexCells[sortperm(cellAngles,rev=true)]
            h = @SVector [0.0,0.0]
            curlSum = 0
            for (i,j) in enumerate(vertexEdges)
                h = h + ϵ*F[k,vertexCells[i]]
                curlSum += A[j,k]*(T[j]⋅h)/linkTriangleAreas[k]
            end
        else
            # # Set angles relative to pressure force angle, equivalent to the angle of a cell that doesn't actually exist
            # pressureAngle = atan((-1.0.*externalF[k])...)
            # vertexEdges = findall(x->x!=0,A[:,k])
            # edgeAngles = zeros(length(vertexEdges))
            # for (i,e) in enumerate(vertexEdges)
            #     edgeAngles[i] = atan((-A[e,k].*edgeTangents[e])...)
            # end
            #
            # edgeAngles .+= 2π-pressureAngle
            # edgeAngles .= edgeAngles.%2π
            # vertexEdges .= vertexEdges[sortperm(edgeAngles,rev=true)]
            #
            # vertexCells = findall(x->x!=0,C[:,k])
            # cellAngles = zeros(length(vertexCells))
            # for i=1:length(cellAngles)
            #     cellAngles[i] = atan((cellPositions[vertexCells[i]].-R[k])...)
            # end
            # cellAngles .+= 2π-pressureAngle
            # cellAngles .= cellAngles.%2π
            # vertexCells .= vertexCells[sortperm(cellAngles,rev=true)]
            #
            # h = @SVector [0.0,0.0]
            # curlSum = 0
            # h = h + ϵ*externalF[k]
            # curlSum += A[vertexEdges[1],k]*(T[vertexEdges[1]]⋅h)/linkTriangleAreas[k]
            # for (i,j) in enumerate(vertexEdges[2:end])
            #     h = h + ϵ*F[k,vertexCells[i]]
            #     curlSum += A[j,k]*(T[j]⋅h)/linkTriangleAreas[k]
            # end
            curlSum = 0.0
        end
        push!(vertexCurls,curlSum)
    end
    return vertexCurls
end

function makeLf(params,matrices,trapeziumAreas)
    @unpack B, Bᵀ, edgeLengths, cellAreas = matrices
    @unpack nCells = params
    onesVec = ones(1,nCells)
    boundaryEdges = abs.(onesVec*B)
    H = Diagonal(cellAreas)
    boundaryEdgesFactor = abs.(boundaryEdges.-1)# =1 for internal vertices, =0 for boundary vertices
    diagonalComponent = (boundaryEdgesFactor'.*((edgeLengths.^2)./(2.0.*trapeziumAreas)))[:,1] # Multiply by boundaryEdgesFactor vector to set boundary vertex contributions to zero
    Tₑ = Diagonal(diagonalComponent)
    Lf = (H\B)*Tₑ*Bᵀ
    dropzeros!(Lf)
    return Lf
end

function makeLc(params,matrices,trapeziumAreas)
    @unpack B, Bᵀ, cellAreas = matrices
    @unpack nCells = params
    onesVec = ones(1,nCells)
    boundaryEdges = abs.(onesVec*B)
    boundaryEdgesFactor = abs.(boundaryEdges.-1)# =1 for internal vertices, =0 for boundary vertices
    H = Diagonal(cellAreas)
    Tₗ = Diagonal(((norm.(T)).^2)./(2.0.*trapeziumAreas))
    invTₗ = inv(Tₗ)
    boundaryEdgesFactorMat = Diagonal(boundaryEdgesFactor[1,:])
    Lc = (H\B)*boundaryEdgesFactorMat*invTₗ*Bᵀ
    dropzeros!(Lc)
    return Lc
end

function makeLv(params,matrices,linkTriangleAreas,trapeziumAreas)
    @unpack A, Aᵀ, edgeLengths = matrices
    E = Diagonal(linkTriangleAreas)
    Tₑ = Diagonal((edgeLengths.^2)./(2.0.*trapeziumAreas))
    Lᵥ = (E\Aᵀ)*(Tₑ\A)
    dropzeros!(Lᵥ)
    return Lᵥ
end

function makeLt(params,matrices,T,linkTriangleAreas,trapeziumAreas)
    @unpack A, Aᵀ = matrices
    E = Diagonal(linkTriangleAreas)
    Tₗ = Diagonal(((norm.(T)).^2)./(2.0.*trapeziumAreas))
    Lₜ = (E\Aᵀ)*Tₗ*A
    dropzeros!(Lₜ)
    return Lₜ
end



function makeCellVerticesDict(params,matrices)
    @unpack C, R, cellPositions = matrices
    @unpack nCells = params
    cellVerticesDict = Dict()
    for i=1:nCells
        # Find vertices around cell
        cellVertices = findall(x->x!=0,C[i,:])
        # Find angles of vertices around cell
        vertexAngles = zeros(size(cellVertices))
        for (k,v) in enumerate(cellVertices)
            vertexAngles[k] = atan((R[v].-cellPositions[i])...)
        end
        # Sort vertices around cell by polar angle
        cellVertices .= cellVertices[sortperm(vertexAngles)]
        # Store sorted cell vertices for this cell
        cellVerticesDict[i] = cellVertices
    end
    return cellVerticesDict
end
