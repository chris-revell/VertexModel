#
#  Visualise.jl
#  VertexModel
#
#  Created by Christopher Revell on 16/02/2021.
#
#
#

module Visualise

# Julia packages
using Printf
using LinearAlgebra
using ColorSchemes
using Colors
using UnPack
using GeometryBasics
using Random
using CairoMakie

function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end

@views function visualise(t,fig,ax,mov,params,matrices)

   plotCells         = 1
   plotEdges         = 0
   scatterEdges      = 1
   scatterVertices   = 1
   scatterCells      = 1

   @unpack R,A,B,Ā,B̄,C,F,cellPositions,edgeTangents,edgeMidpoints,boundaryVertices,vertexEdges = matrices
   @unpack nEdges,nVerts,nCells = params

   empty!(ax)

   ax.title="t = $t"

   # Plot cells
   if plotCells == 1
      for i=1:nCells
         cellVertices = findall(x->x!=0,C[i,:])
         vertexAngles = zeros(size(cellVertices))
         for (k,v) in enumerate(cellVertices)
            vertexAngles[k] = atan((R[v].-cellPositions[i])...)
         end
         cellVertices .= cellVertices[sortperm(vertexAngles)]
         poly!(ax,Point2f.(R[cellVertices]),color=(getRandomColor(i),0.5))
      end
   end

   # Scatter vertices
   if scatterEdges == 1
      scatter!(ax,Point2f.(R),color=:green)
      annotations!(ax,string.(collect(1:length(R))), Point2f.(R),color=:green)
   end

   # Scatter edge midpoints
   if scatterEdges == 1
      scatter!(ax,Point2f.(edgeMidpoints),color=:blue)
      annotations!(ax,string.(collect(1:length(edgeMidpoints))), Point2f.(edgeMidpoints),color=:blue)
   end

   # Scatter cell positions
   if scatterCells == 1
      scatter!(ax,Point2f.(cellPositions),color=:red)
      annotations!(ax,string.(collect(1:length(cellPositions))), Point2f.(cellPositions),color=:red)
   end

   # Plot edges
   # For each edge, use Ā adjacency matrix to find corresponding vertices x, and plot line between x[1] and x[2]
   if plotEdges == 1
      xs = Point2f[]
      us = Vec2f[]
      colours = Tuple{Symbol, Float64}[]
      for c=1:nCells
         es = findall(x->x!=0,B[c,:])
         for i in es
            vs =findall(x->x!=0,Ā[i,:])
            colour=:black
            # Use B to set colour of edge depending on whether it runs with or against the orientation of the cell face
            if matrices.B[c,i] < 0
               colour = (:red,0.25)
            else
               colour = (:blue,0.25)
            end
            # Use A to set direction of arrow along edge
            if A[i,vs[1]] < 0
               #arrows!(ax,Point2f.([R[x[1]]]), Vec2f.([edgeTangents[i]]),color=(colour,0.25),arrowsize=25,linewidth=5)
               push!(xs, Point2f(R[vs[1]]))
               push!(us, Vec2f(edgeTangents[i]))
               push!(colours, colour)
            else
               #arrows!(ax,Point2f.([R[x[2]]]), Vec2f.([edgeTangents[i]]),color=(colour,0.25),arrowsize=25,linewidth=5)
               push!(xs, Point2f(R[vs[2]]))
               push!(us, Vec2f(edgeTangents[i]))
               push!(colours, colour)
            end
         end
      end
      arrows!(ax,xs,us,color=colours,arrowsize=25,linewidth=5)
   end

   recordframe!(mov)

   return 0

end

export visualise

end
