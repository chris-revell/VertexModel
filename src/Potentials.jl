#
#  Potentials.jl
#  VertexModel
#
#  Created by Christopher Revell on 14/02/2021.
#
#

module Potentials

# Julia packages
using DrWatson
using LinearAlgebra
using SparseArrays
using FromFile
using GeometryBasics

# Local modules
@from "$(projectdir("src","AnalysisFunctions.jl"))" using AnalysisFunctions
@from "$(projectdir("src","Laplacians.jl"))" using Laplacians

function psicPotential(R,params,matrices)
    T = makeCellLinks(params,matrices)
    edgeTrapezia = makeEdgeTrapezia(R,params,matrices)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    linkTriangles = makeLinkTriangles(R,params,matrices)
    linkTriangleAreas = abs.(area.(linkTriangles))
    Lf = makeLf(params,matrices,trapeziumAreas)
    cellDivs = -1.0.*calculateCellDivs(R,params,matrices)
    onesVec = ones(params.nCells)
    H = Diagonal(matrices.cellAreas)
    eigenvectors = (eigen(Matrix(Lf))).vectors
    eigenvalues = (eigen(Matrix(Lf))).values
    ḡ = ((onesVec'*H*cellDivs)/(onesVec'*H*ones(params.nCells))).*onesVec
    ğ = cellDivs.-ḡ
    ψ̆ = zeros(params.nCells)
    spectrum = Float64[]
    for k=2:params.nCells
        numerator = eigenvectors[:,k]'*H*ğ
        denominator = eigenvalues[k]*(eigenvectors[:,k]'*H*eigenvectors[:,k])
        ψ̆ .-= (numerator/denominator).*eigenvectors[:,k]
        push!(spectrum,(numerator/denominator))
    end
    return ψ̆, spectrum
end

function psivPotential(R,params,matrices)
    T = makeCellLinks(params,matrices)
    edgeTrapezia = makeEdgeTrapezia(R,params,matrices)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    linkTriangles = makeLinkTriangles(R,params,matrices)
    linkTriangleAreas = abs.(area.(linkTriangles))
    q = calculateSpokes(R,params,matrices)
    Lₜ = makeLt(params,matrices,T,linkTriangleAreas,trapeziumAreas)
    eigenvectors = (eigen(Matrix(Lₜ))).vectors
    eigenvalues = (eigen(Matrix(Lₜ))).values
    vertexDivs = -1.0.*calculateVertexDivs(R,params,matrices,q,linkTriangleAreas)
    onesVec = ones(params.nVerts)
    E = Diagonal(linkTriangleAreas)
    ḡ = ((onesVec'*E*vertexDivs)/(onesVec'*E*ones(params.nVerts))).*onesVec
    ğ = vertexDivs.-ḡ
    ψ̆ = zeros(params.nVerts)
    spectrum = Float64[]
    for k=2:params.nVerts
        numerator = -eigenvectors[:,k]'*E*ğ
        denominator = eigenvalues[k]*(eigenvectors[:,k]'*E*eigenvectors[:,k])
        ψ̆ .+= (numerator/denominator).*eigenvectors[:,k]
        push!(spectrum,(numerator/denominator))
    end
    return ψ̆, spectrum
end

function capitalPsivPotential(R,params,matrices)
    T = makeCellLinks(params,matrices)
    edgeTrapezia = makeEdgeTrapezia(R,params,matrices)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    linkTriangles = makeLinkTriangles(R,params,matrices)
    linkTriangleAreas = abs.(area.(linkTriangles))
    Lₜ = makeLt(params,matrices,T,linkTriangleAreas,trapeziumAreas)
    eigenvectors = (eigen(Matrix(Lₜ))).vectors
    eigenvalues = (eigen(Matrix(Lₜ))).values
    q = calculateSpokes(R,params,matrices)
    vertexCurls = calculateVertexCurls(R,params,matrices,q,linkTriangleAreas)
    onesVec = ones(params.nVerts)
    E = Diagonal(linkTriangleAreas)
    ḡ = ((onesVec'*E*vertexCurls)/(onesVec'*E*ones(params.nVerts))).*onesVec
    ğ = vertexCurls.-ḡ
    ψ̆ = zeros(params.nVerts)
    spectrum = Float64[]
    for k=2:params.nVerts
        numerator = -eigenvectors[:,k]'*E*ğ
        denominator = eigenvalues[k]*(eigenvectors[:,k]'*E*eigenvectors[:,k])
        ψ̆ .+= (numerator/denominator).*eigenvectors[:,k]
        push!(spectrum,(numerator/denominator))
    end
    return ψ̆, spectrum
end

export psicPotential, psivPotential, capitalPsivPotential

end #end module 