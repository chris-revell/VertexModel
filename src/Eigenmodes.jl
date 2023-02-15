#
#  Eigenmodes.jl
#  VertexModel
#
#  Created by Christopher Revell on 14/02/2021.
#
#

module Eigenmodes

# Julia packages
using LinearAlgebra
using SparseArrays
using FromFile
using GeometryBasics

# Local modules
@from "AnalysisFunctions" using AnalysisFunctions
@from "Laplacians" using Laplacians

function eigenmodesLt(matrices,params)    
    T = makeCellLinks(params,matrices)
    edgeTrapezia = makeEdgeTrapezia(params,matrices)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    linkTriangles = makeLinkTriangles(params,matrices)
    linkTriangleAreas = abs.(area.(linkTriangles))
    Lₜ = makeLt(params,matrices,T,linkTriangleAreas,trapeziumAreas)
    decomposition = (eigen(Matrix(Lₜ))).vectors
    return decomposition
end

function allEigenModesLf(matrices,params)
    T = makeCellLinks(params,matrices)
    edgeTrapezia = makeEdgeTrapezia(params,matrices)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    linkTriangles = makeLinkTriangles(params,matrices)
    linkTriangleAreas = abs.(area.(linkTriangles))
    Lf = makeLf(params,matrices,trapeziumAreas)
    decomposition = (eigen(Matrix(Lf))).vectors
    return decomposition
end

function eigenmodelLv(matrices,params)
    T = makeCellLinks(params,matrices)
    edgeTrapezia = makeEdgeTrapezia(params,matrices)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    linkTriangles = makeLinkTriangles(params,matrices)
    linkTriangleAreas = abs.(area.(linkTriangles))    
    Lᵥ = makeLv(params,matrices,linkTriangleAreas,trapeziumAreas)
    decomposition = (eigen(Matrix(Lᵥ))).vectors
end 

export eigenmodesLt, eigenmodesLf, eigenmodesLv

end #end module 