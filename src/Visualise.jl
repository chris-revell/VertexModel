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

@views function visualise(fig,ax,mov,params,matrices)


   @unpack R,A,B,Ā,B̄,C,F,cellPositions,edgeTangents,edgeMidpoints,boundaryVertices,vertexEdges = matrices
   @unpack nEdges,nVerts,nCells = params

   empty!(ax)

   plotCells         = 1
   plotEdges         = 1
   scatterEdges      = 1
   scatterVertices   = 1
   scatterCells      = 1

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
            x=findall(x->x!=0,Ā[i,:])
            colour=:black
            # Use B to set colour of edge depending on whether it runs with or against the orientation of the cell face
            if matrices.B[c,i] < 0
               colour = (:red,0.25)
            else
               colour = (:blue,0.25)
            end
            # Use A to set direction of arrow along edge
            if A[i,x[1]] < 0
               #arrows!(ax,Point2f.([R[x[1]]]), Vec2f.([edgeTangents[i]]),color=(colour,0.25),arrowsize=25,linewidth=5)
               push!(xs, Point2f(R[x[1]]))
               push!(us, Vec2f(edgeTangents[i]))
               push!(colours, colour)
            else
               #arrows!(ax,Point2f.([R[x[2]]]), Vec2f.([edgeTangents[i]]),color=(colour,0.25),arrowsize=25,linewidth=5)
               push!(xs, Point2f(R[x[2]]))
               push!(us, Vec2f(edgeTangents[i]))
               push!(colours, colour)
            end
         end
      end
      arrows!(ax,xs,us,color=colours,arrowsize=25,linewidth=5)
   end

   recordframe!(mov)

   # Vertex moment kites
   # for i=1:nVerts
   #    # Exclude boundary vertices.
   #    if boundaryVertices[i] == 1
   #       # Do nothing
   #    else
   #       vertexEdges = findall(x->x!=0,A[:,i])
   #       # Loop over 3 edges around the vertex.
   #       for k=0:2
   #          # Find vector separating the midpoint of an edge and the midpoint of the next edge around the vertex.
   #          dM = edgeMidpoints[vertexEdges[i,(k+1)%3+1]] .- edgeMidpoints[vertexEdges[i,k+1]] # Could use arrayLoop function here
   #          # Find cell bordered by both these two edges
   #          cellID = intersect(findall(x->x!=0,B̄[:,vertexEdges[i,(k+1)%3+1]]),findall(x->x!=0,B̄[:,vertexEdges[i,k+1]]))[1]
   #          kite = Shape(Point2f.([cellPositions[cellID],edgeMidpoints[vertexEdges[i,k+1]],R[i],edgeMidpoints[vertexEdges[i,(k+1)%3+1]]]))
   #          dotProduct = normalize(edgeTangents[vertexEdges[i,((k+2)%3)+1]])⋅normalize(dM)
   #          plot!(kite,linewidth=0,fillcolor=:blue,fillalpha=abs(dotProduct),seriestype=:shape)
   #       end
   #    end
   # end

   return 0

end

export visualise

end
