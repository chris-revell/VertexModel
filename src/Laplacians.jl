#
#  Laplacians.jl
#  VertexModel
#
#  Created by Christopher Revell on 15/02/2021.
#
#

module Laplacians

# Julia packages
using DrWatson
using LinearAlgebra
using SparseArrays
using FromFile

# Local modules
@from "AnalysisFunctions.jl" using AnalysisFunctions

function makeLf(params,matrices,trapeziumAreas)
    @unpack B,Bᵀ,cellAreas,edgeLengths = matrices
    @unpack nCells = params
    onesVec = ones(1,nCells)
    boundaryEdges = abs.(onesVec*B)
    H = Diagonal(cellAreas)
    boundaryEdgesFactor = abs.(boundaryEdges.-1)# =1 for internal vertices, =0 for boundary vertices
    diagonalComponent = (boundaryEdgesFactor'.*((edgeLengths.^2)./(2.0.*trapeziumAreas)))[:,1] # Multiply by boundaryEdgesFactor vector to set boundary vertex contributions to zero
    Tₑ = Diagonal(diagonalComponent)
    invH = inv(H)
    Lf = invH*B*Tₑ*Bᵀ
    dropzeros!(Lf)    
    return Lf
end

function makeLc(params,matrices,T,trapeziumAreas)
    @unpack B,Bᵀ,cellAreas = matrices
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
    @unpack A,Aᵀ,edgeLengths = matrices
    E = Diagonal(linkTriangleAreas)
    Tₑ = Diagonal((edgeLengths.^2)./(2.0.*trapeziumAreas))
    Lᵥ = (E\Aᵀ)*(Tₑ\A)
    dropzeros!(Lᵥ)
    return Lᵥ
end

function makeLt(params,matrices,T,linkTriangleAreas,trapeziumAreas)
    @unpack A,Aᵀ = matrices
    E = Diagonal(linkTriangleAreas)
    Tₗ = Diagonal(((norm.(T)).^2)./(2.0.*trapeziumAreas))
    Lₜ = (E\Aᵀ)*Tₗ*A
    dropzeros!(Lₜ)
    return Lₜ
end

export makeLf, makeLc, makeLv, makeLt

end #end module 