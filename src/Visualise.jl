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
using Plots
using Printf
using LinearAlgebra
using ColorSchemes
using UnPack
using GeometryBasics
using Random
#using CairoMakie

function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end

@views function visualise(anim,params,matrices)


   @unpack R,A,B,Ā,B̄,C,F,cellPositions,edgeTangents,edgeMidpoints,boundaryVertices,vertexEdges = matrices
   @unpack nEdges,nVerts,nCells = params

   # Create plot canvas
   plot(aspect_ratio=:equal,border=:none,legend=:false,dpi=300,size=(500,500))

   # Plot cells
   for i=1:nCells
      cellVertices = findall(x->x!=0,C[i,:])
      vertexAngles = zeros(size(cellVertices))
      for (k,v) in enumerate(cellVertices)
         vertexAngles[k] = atan((R[v].-cellPositions[i])...)
      end
      cellVertices .= cellVertices[sortperm(vertexAngles)]
      Cell = Shape(Point2f.(R[cellVertices]))
      plot!(Cell,color=getRandomColor(i),alpha=0.5,linealpha=0.0)
   end

   # Scatter vertices
   # scatter!(Point2f.(R),series_annotations=text.(1:length(R),:bottom),markersize=2,seriescolor=:green)

   # # Scatter edge midpoints
   # scatter!(Point2f.(edgeMidpoints),color=:white,markersize=1,series_annotations=text.(1:length(edgeMidpoints),:bottom),seriescolor=:blue)

   # Scatter cell positions
   # scatter!(Point2f.(cellPositions),color=:red,markersize=1,markerstroke=:red,series_annotations=text.(1:length(cellPositions),:bottom),seriescolor=:red)

   # Plot edges
   # For each edge, use Ā adjacency matrix to find corresponding vertices x, and plot line between x[1] and x[2]
   # for c=1:nCells
   #    es = findall(x->x!=0,B[c,:])
   #    for i in es
   #       x=findall(x->x!=0,Ā[i,:])
   #       colour=:black
   #       # Use B to set colour of edge depending on whether it runs with or against the orientation of the cell face
   #       if matrices.B[c,i] < 0
   #          colour = :red
   #       else
   #          colour = :blue
   #       end
   #       # Use A to set direction of arrow along edge
   #       if A[i,x[1]] < 0
   #          plot!(Point2f.([R[x[1]], R[x[2]]]),arrow=true,color=colour,linewidth=2,arrowsize=16,alpha=0.25)
   #       else
   #          plot!(Point2f.([R[x[2]], R[x[1]]]),arrow=true,color=colour,linewidth=2,arrowsize=16,alpha=0.25)
   #       end
   #    end
   # end

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

   frame(anim)

   return 0

end

export visualise

end
