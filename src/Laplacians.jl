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
using StaticArrays
using FromFile
using BlockArrays

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

function makeG(params)
    @unpack nCells, γ = params
    G=Diagonal(vcat(fill(1,nCells), fill(γ, nCells)))
    return G
end

function makeM(matrices)
    @unpack A, Ā, B, B̄, edgeTangents, ϵ = matrices
    #dAdr= - 1/2 Sum j ϵᵢ . Bᵢⱼ Āⱼₖ tⱼ = 1/2 Sum j Ānᵢⱼ = 1/2 B diag(ϵ.t) Ā
    dAdr=-1/2*(B*Diagonal(eachcol((ϵ*reduce(vcat,transpose(edgeTangents))')))*Ā)
    #dLdr=Sum j B̄ᵢⱼ Aⱼₖ t̂ⱼ   = B̄ diag(t̂) A
    dLdr= B̄*Diagonal(edgeTangents./norm.(edgeTangents))*A
    
    M=vcat(dAdr,dLdr)
    dropzeros!(M)
  
    return M
end


function makeEvLc(M)
    EvLc=tr.(M*M')
    dropzeros!(EvLc)
    return EvLc
end

function makeEvLv(M,G)
    EvLv=sparse([x' for x in M']*[x' for x in (G*M)])
    dropzeros!(EvLv)
    return EvLv
end

function makeX(params,matrices)
   
    @unpack A, Ā, B, B̄, edgeLengths, edgeTangents, ϵ = matrices
    @unpack nCells, nVerts = params
   
    dAdrr = Array{SMatrix{2,2,Float64}}(undef,nCells,nVerts,nVerts)
    fill!(dAdrr,@SMatrix zeros(2,2))
    dLdrr = Array{SMatrix{2,2,Float64}}(undef,nCells,nVerts,nVerts)
    fill!(dLdrr,@SMatrix zeros(2,2))

    n=eachcol((ϵ*reduce(vcat,transpose(edgeTangents))'))

    for k=1:nVerts 
        for m=1:nVerts
            for j in nzrange(A,k)
                for i in nzrange(B,rowvals(A)[j])

                    dAdrr[rowvals(B)[i], k,m ]+= -0.5*B[rowvals(B)[i],rowvals(A)[j]]*Ā[rowvals(A)[j],k]*A[rowvals(A)[j],m].*(ϵ)

                    dLdrr[rowvals(B)[i], k,m ]+=(B̄[rowvals(B)[i],rowvals(A)[j]]*A[rowvals(A)[j],k]*A[rowvals(A)[j],m]*(n[rowvals(A)[j]]/edgeLengths[rowvals(A)[j]])*(n[rowvals(A)[j]]/edgeLengths[rowvals(A)[j]])')/edgeLengths[rowvals(A)[j]]

                end

            end
        end
    end

    X=vcat(dAdrr, dLdrr)

    return X
end



function makeD(params,matrices,X, Lvevals, Lvevec, q)
    #assuming LAPACK.syev used to find eigenvectors of Lv
    #q is number of non zero eigen vectors, generally 2*nCells -1

    #Cell tensions are negative due to previous sign conventions in code.
    
    @unpack cellTensions, cellPressures = matrices
    @unpack nCells, nVerts = params
    g=vcat(cellPressures, -cellTensions)

    gX=Matrix{SMatrix{2,2,Float64,4}}(undef,nVerts,nVerts)
    fill!(gX,@SMatrix zeros(2,2))
    for α=1:2*nCells
        gX+=g[α]X[α, :,:]
    end

    D=zeros(q,q)
    qeval=Lvevals[2*nVerts+1-q: 2*nVerts]
    qevec=Lvevec[:, 2*nVerts+1-q:2*nVerts]
    
    D=qevec'*Matrix(mortar(gX))*qevec+ Diagonal(qeval)

    return D

end

export makeLf, makeLc, makeLv, makeLt, makeG, makeM, makeEvLc, makeEvLv, makeX, makeD

end #end module 