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
   #plot(xlims=(-0.5,0.5),ylims=(-0.5,0.5),aspect_ratio=:equal,color=:black,legend=:false,border=:none,markersize=4,markerstroke=:black,dpi=300,size=(1000,1000))
   plot(aspect_ratio=:equal,color=:black,legend=:false,border=:none,markersize=4,markerstroke=:black,dpi=300,size=(1000,1000))

   # Scatter vertices
   scatter!(Point2f.(R),series_annotations=text.(1:length(R),:bottom),markersize=10)

   # # Scatter edge midpoints
   scatter!(Point2f.(edgeMidpoints),color=:white,markersize=4,series_annotations=text.(1:length(edgeMidpoints),:bottom))

   # Scatter cell positions
   scatter!(Point2f.(cellPositions),color=:red,markersize=4,markerstroke=:red,series_annotations=text.(1:length(cellPositions),:bottom))

   # Plot edges
   # For each edge, use Ā adjacency matrix to find corresponding vertices x, and plot line between x[1] and x[2]
   for c=1:nCells
      es = findall(x->x!=0,B[c,:])
      for i in es
         x=findall(x->x!=0,Ā[i,:])
         colour=:black
         if matrices.B[c,i] < 0
            colour = :red
         else
            colour = :blue
         end
         if A[i,x[1]] < 0
            #arrows!([R[x[1]]],[R[x[2]]],[edgeTangents[i][1]],[edgeTangents[i][2]],color=colour)
            plot!(Point2f.([R[x[1]], R[x[2]]]),arrow=true,color=colour,linewidth=8,arrowsize=16,alpha=0.5)
         else
            #arrows!([R[x[1]]],[R[x[2]]],[edgeTangents[i][1]],[edgeTangents[i][2]],color=colour)
            plot!(Point2f.([R[x[2]], R[x[1]]]),arrow=true,color=colour,linewidth=8,arrowsize=16,alpha=0.5)
         end
      end
   end

   # Plot cells
   for i=1:nCells
      cellVertices = findall(x->x!=0,C[i,:])
      vertexAngles = zeros(size(cellVertices))
      for (k,v) in enumerate(cellVertices)
         vertexAngles[k] = atan((R[v].-cellPositions[i])...)
      end
      cellVertices .= cellVertices[sortperm(vertexAngles)]
      Cell = Shape(Point2f.(R[cellVertices]))
      plot!(Cell,color=getRandomColor(i),alpha=0.5)
   end

   # Vertex moment kites
   # for i=1:nVerts
   #    # Exclude boundary vertices.
   #    if boundaryVertices[i] == 0
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
   #    else
   #       # Skip boundary vertices
   #    end
   # end


   # *************** Force vectors ***************************
   # xs = zeros(0)
   # ys = zeros(0)
   # us = zeros(0)
   # vs = zeros(0)
   # for i=1:nVerts
   #    x=findall(x->x>0.5,C[:,i])
   #    for j in x
   #       append!(xs,R[i,1])
   #       append!(ys,R[i,2])
   #       append!(us,5.0*F[i,j,2])
   #       append!(vs,-5.0*F[i,j,1])
   #    end
   # end
   # quiver!(xs,ys,quiver=(us,vs),color=:orange,legend=:false)
   # *************** End force vectors ***********************

   frame(anim)

end

export visualise

end
