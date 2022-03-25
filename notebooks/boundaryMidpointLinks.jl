# Import Julia packages
using DrWatson
@quickactivate
using Revise
using LinearAlgebra
using DelimitedFiles
using SparseArrays
using StaticArrays
using CairoMakie
using UnPack
using GeometryBasics
using Random
using Colors
using JLD2

# Local modules
includet("$(projectdir())/src/VertexModelContainers.jl"); using .VertexModelContainers
includet("$(projectdir())/notebooks/functions.jl")

#dataDirectory = "data/sims/2022-02-28-19-30-22"

# Import system data
conditionsDict    = load("$dataDirectory/dataFinal.jld2")
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$dataDirectory/matricesFinal.jld2")
@unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths,externalF = matricesDict["matrices"]

onesVec = ones(1,nCells)
boundaryEdges = abs.(onesVec*B)
cᵖ = boundaryEdges'.*edgeMidpoints

s = Matrix{SVector{2,Float64}}(undef,nCells,nVerts)
p = Matrix{SVector{2,Float64}}(undef,nCells,nVerts)
fill!(s,@SVector zeros(2))
fill!(p,@SVector zeros(2))
for i=1:nCells
    for k=1:nVerts
        for j=1:nEdges
            s[i,k] += B[i,j].*cᵖ[j].*A[j,k]
        end
        p[i,k] = ϵ*s[i,k]
    end
end

s_reduced = s[findall(!iszero,s)]

for k=1:nVerts
    if boundaryVertices[k] != 0
        vertexEdges = findall(!iszero,A[:,k])
        if length(vertexEdges) == 2
            rotatedExternalForce = ϵ*externalF[k]
            edgeCell = findall(!iszero,C[:,k])
            println("s")
            display(2.0*pressureExternal.*s[edgeCell[1],k])
            println("rotated external force ")
            display(rotatedExternalForce)
        end
    end
end

function calculateVertexDivs(params,matrices,T,linkTriangleAreas)
    @unpack A, C, R, edgeTangents, ϵ, boundaryVertices = matrices
    @unpack nVerts = params

    onesVec = ones(1,nCells)
    boundaryEdges = abs.(onesVec*B)
    cᵖ = boundaryEdges'.*edgeMidpoints

    s = Matrix{SVector{2,Float64}}(undef,nCells,nVerts)
    p = Matrix{SVector{2,Float64}}(undef,nCells,nVerts)
    fill!(s,@SVector zeros(2))
    fill!(p,@SVector zeros(2))
    for i=1:nCells
        for k=1:nVerts
            for j=1:nEdges
                s[i,k] += B[i,j].*cᵖ[j].*A[j,k]
            end
            p[i,k] = ϵ*s[i,k]
        end
    end

    # Rotation matrix around vertices is the opposite of that around cells
    ϵₖ = -1*ϵ

    vertexDivs = Float64[]
    # Working around a given vertex, an h force space point from a cell is mapped to the next edge anticlockwise from the cell
    for k=1:nVerts
        if boundaryVertices[k] == 0
            # Internal vertices

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
            # Peripheral vertices
            vertexCells = findall(x->x!=0,C[:,k])

            if length(vertexCells)==1
                # Peripheral vertex with a single kite

                cellAngle = atan((cellPositions[vertexCells[1]].-R[k])...)

                vertexEdges = findall(x->x!=0,A[:,k])
                edgeAngles = zeros(length(vertexEdges))
                for (i,e) in enumerate(vertexEdges)
                    edgeAngles[i] = atan((-A[e,k].*edgeTangents[e])...)
                end
                edgeAngles .= edgeAngles .+ 2π .- cellAngle
                edgeAngles .= edgeAngles.%2π
                vertexEdges .= vertexEdges[sortperm(edgeAngles,rev=true)]

                splitExternalF = SVector{2,Float64}[]
                push!(splitExternalF, 0.5*pressureExternal*B[vertexCells[1],vertexEdges[1]]*(ϵ*edgeTangents[vertexEdges[1]]))
                push!(splitExternalF, 0.5*pressureExternal*B[vertexCells[1],vertexEdges[2]]*(ϵ*edgeTangents[vertexEdges[2]]))

                h = @SVector [0.0,0.0]
                divSum = 0

                hⱼ = h + ϵ*F[k,vertexCells[1]]
                divSum -= A[vertexEdges[1],k]*((ϵₖ*T[vertexEdges[1]])⋅hⱼ)/linkTriangleAreas[k]
                hₖ = hⱼ .+ ϵ*splitExternalF[1]
                divSum += p[i,k]⋅hₖ


                for (i,j) in enumerate(vertexEdges[2:end])
                    h = h + ϵ*F[k,vertexCells[i]]
                    divSum -= A[j,k]*((ϵₖ*T[j])⋅h)/linkTriangleAreas[k]
                end


            elseif length(vertexCells)==2
                # Peripheral vertex with 2 kites

                cellAngles = zeros(length(vertexCells))
                for i=1:length(cellAngles)
                    cellAngles[i] = atan((cellPositions[vertexCells[i]].-R[k])...)
                end

                vertexEdges = findall(x->x!=0,A[:,k])
                edgeAngles = zeros(length(vertexEdges))
                for (i,e) in enumerate(vertexEdges)
                    edgeAngles[i] = atan((-A[e,k].*edgeTangents[e])...)
                end


                edgeAngles .= edgeAngles .+ 2π .- cellAngle
                edgeAngles .= edgeAngles.%2π
                vertexEdges .= vertexEdges[sortperm(edgeAngles,rev=true)]

                splitExternalF = SVector{2,Float64}[]
                push!(splitExternalF, 0.5*pressureExternal*B[vertexCells[1],vertexEdges[1]]*(ϵ*edgeTangents[vertexEdges[1]]))
                push!(splitExternalF, 0.5*pressureExternal*B[vertexCells[1],vertexEdges[2]]*(ϵ*edgeTangents[vertexEdges[2]]))

                h = @SVector [0.0,0.0]
                divSum = 0

                h = h + ϵ*F[k,vertexCells[1]]
                divSum -= A[vertexEdges[1],k]*((ϵₖ*T[vertexEdges[1]])⋅h)/linkTriangleAreas[k]



                for (i,j) in enumerate(vertexEdges[2:end])
                    h = h + ϵ*F[k,vertexCells[i]]
                    divSum -= A[j,k]*((ϵₖ*T[j])⋅h)/linkTriangleAreas[k]
                end






        end
        push!(vertexDivs,divSum)
    end
    return vertexDivs
end
