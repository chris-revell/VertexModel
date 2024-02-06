#
#  CellProperties.jl
#  VertexModel
#
#  Created by Natasha Cowley on 14/11/2023.
#
#

module CellProperties

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


function makeShapeTensors(R,params,matrices)
    @unpack C, cellEdgeCount, cellPositions= matrices
    @unpack nCells = params
    
    cellShapeTensors= Array{SMatrix{2,2,Float64}}(undef,nCells)
    fill!(cellShapeTensors,@SMatrix zeros(2,2))

    for c=1:nCells
        cellVertices = findall(x->x!=0,C[c,:])
        Rα=R[cellVertices].-cellPositions[c, :]
        cellShapeTensors[c]=tr(Rα*Rα')./cellEdgeCount[c]
    end

    return cellShapeTensors
end

function getPeff(params, matrices)
    @unpack cellAreas, cellPerimeters, cellPressures, cellTensions= matrices
    @unpack nCells, γ, L₀= params

    Peff=zeros(nCells)
    Peff=cellPressures+(-cellTensions.*cellPerimeters)./(2*cellAreas)
    return Peff
end


function makeCellQandJ(params, matrices)
    @unpack B, edgeTangents, edgeLengths,  cellPerimeters= matrices
    @unpack nCells= params

    cellQ= fill(SMatrix{2,2,Float64}(zeros(2,2)), nCells)
    cellJ= fill(SMatrix{2,2,Float64}(zeros(2,2)), nCells)

    


    for c=1:nCells
        cellEdges = findall(x->x!=0,B[c,:])
        cellUnitTangents=edgeTangents[cellEdges]./edgeLengths[cellEdges]
        cellQ[c]=tr(edgeTangents[cellEdges]*cellUnitTangents')./cellPerimeters[c]
        cellJ[c]=cellQ[c]-0.5*I
    end

    return cellQ, cellJ
end

function getShearStress(params, matrices, cellJ)
    @unpack cellAreas, cellTensions, cellPerimeters= matrices
    @unpack nCells= params

    detJ=det.(cellJ)
    cellShearStress=((cellPerimeters.*(-cellTensions))./cellAreas).*sqrt.(-detJ)

    return cellShearStress
end



function getCircularity(params, cellShapeTensors)
@unpack nCells = params
cellCircularity=zeros(nCells)
for c=1:nCells
    λ₁, λ₂= eigen(cellShapeTensors[c]).values
    cellCircularity[c]=abs( λ₁ / λ₂)
end

return cellCircularity

end

function getShapeStressAngle(params, matrices, cellShapeTensors, Peff)
    #from 0 to pi
    @unpack nCells= params

    shapeAngle=zeros(nCells)
    stressAngle=zeros(ncells)

    for c=1:nCells
        eigenvectors=Matrix(qr(eigen(cellShapeTensors[c]).vectors).Q)
        v₁ = eigenvectors[:,1]
        v₂ = eigenvectors[:,2]
        shapeAngle[c]=atan(v₂[2], v₂[1])
        if shapeAngle[c]<0; shapeAngle[c] +=pi end

        if Peff[c]<0
            stressAngle[c]=atan(v₁[2], v₁[1])
            if StressAngle[c]<0; ShapeAngle[c] +=pi end
        else
            stressAngle[c]=shapeAngle[c]
        end

    end

    return shapeAngle, stressAngle
end



export makeShapeTensors, getPeff, makeCellQandJ, getCircularity, getShapeStressAngle, getShearStress

end