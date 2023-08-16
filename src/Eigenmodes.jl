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
using DrWatson

# Local modules
@from "$(srcdir("AnalysisFunctions.jl"))" using AnalysisFunctions
@from "$(srcdir("Laplacians.jl"))" using Laplacians

function eigenmodesLt(R,matrices,params)    
    T = makeCellLinks(params,matrices)
    edgeTrapezia = makeEdgeTrapezia(R,params,matrices)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    linkTriangles = makeLinkTriangles(R,params,matrices)
    linkTriangleAreas = abs.(area.(linkTriangles))
    Lₜ = makeLt(params,matrices,T,linkTriangleAreas,trapeziumAreas)
    decomposition = (eigen(Matrix(Lₜ))).vectors
    return decomposition
end

function eigenmodesLf(R,matrices,params)
    edgeTrapezia = makeEdgeTrapezia(R,params,matrices)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    Lf = makeLf(params,matrices,trapeziumAreas)
    decomposition = (eigen(Matrix(Lf))).vectors
    return decomposition
end

function eigenmodesLv(R,matrices,params)
    edgeTrapezia = makeEdgeTrapezia(R,params,matrices)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    linkTriangles = makeLinkTriangles(R,params,matrices)
    linkTriangleAreas = abs.(area.(linkTriangles))    
    Lᵥ = makeLv(params,matrices,linkTriangleAreas,trapeziumAreas)
    decomposition = (eigen(Matrix(Lᵥ))).vectors
end 

export eigenmodesLt, eigenmodesLf, eigenmodesLv

end #end module 