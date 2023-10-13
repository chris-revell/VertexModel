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


function makeM(params,matrices)
    @unpack A, Ā, Aᵀ, B, B̄, Bᵀ,cellAreas,edgeLengths, edgeTangents= matrices
    @unpack nCells, nEdges, nVerts = params
    #dAdr= - 1/2 Sum j ϵᵢ . Bᵢⱼ Āⱼₖ tⱼ = 1/2 Sum j Ānᵢⱼ = 1/2 B diag(ϵ.t) Ā
    dAdr=-1/2*(B*Diagonal(eachcol((ϵ*reduce(vcat,transpose(t))')))*Ā)
    #dLdr=Sum j B̄ᵢⱼ Aⱼₖ t̂ⱼ   = B̄ diag(t̂) A
    dLdr= B̄*Diagonal(t./norm.(t))*A
    
    M=vcat(dAdr,dLdr)
  
    return M
end


function makeL(M)
    L=tr.(M*M')
    return L


export makeLf, makeLc, makeLv, makeLt

end #end module 