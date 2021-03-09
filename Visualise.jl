#
#  Visualise.jl
#  VertexModelJL
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

# Local modules

@inline @views function visualise(A,Ā,B̄,R,C,F,cellPositions,edgeTangents,edgeMidpoints,nEdges,nVerts,nCells,outputCount,folderName,ϵ,boundaryVertices,vertexEdges)

   # Create plot canvas
   plot(xlims=(-2,2),ylims=(-2,2),aspect_ratio=:equal,color=:black,legend=:false,border=:none,markersize=4,markerstroke=:black,dpi=300,size=(1000,1000))

   # Scatter vertices
   scatter!(R[:,1],R[:,2])

   # # Scatter edge midpoints
   # scatter!(edgeMidpoints[:,1],edgeMidpoints[:,2],color=:white,markersize=4)

   # Scatter cell positions
   scatter!(cellPositions[:,1],cellPositions[:,2],color=:red,markersize=4,markerstroke=:red)

   # Plot edges
   # For each edge, use Ā adjacency matrix to find corresponding vertices x, and plot line between x[1] and x[2]
   for i=1:nEdges
      x=findall(x->x!=0,Ā[i,:])
      plot!([R[x[1],1],R[x[2],1]],[R[x[1],2],R[x[2],2]],color=:black,linewidth=4)
   end

   # Vertex moment kites
   for i=1:nVerts
      # Exclude boundary vertices.
      if boundaryVertices[i] == 0
         # Loop over 3 edges around the vertex.
         for k=0:2
            # Find vector separating the midpoint of an edge and the midpoint of the next edge around the vertex.
            dM = edgeMidpoints[vertexEdges[i,(k+1)%3+1],:] .- edgeMidpoints[vertexEdges[i,k+1],:]
            # Find cell bordered by these two edges
            cellID = intersect(findall(x->x!=0,B̄[:,vertexEdges[i,(k+1)%3+1]]),findall(x->x!=0,B̄[:,vertexEdges[i,k+1]]))[1]
            xs = [cellPositions[cellID,1],edgeMidpoints[vertexEdges[i,k+1],1],R[i,1],edgeMidpoints[vertexEdges[i,(k+1)%3+1],1],cellPositions[cellID,1]]
            ys = [cellPositions[cellID,2],edgeMidpoints[vertexEdges[i,k+1],2],R[i,2],edgeMidpoints[vertexEdges[i,(k+1)%3+1],2],cellPositions[cellID,2]]
            dotProduct = normalize(edgeTangents[vertexEdges[i,((k+2)%3)+1],:])⋅normalize(dM)
            plot!(xs,ys,linewidth=0,fillcolor=:blue,fillalpha=abs(dotProduct),seriestype=:shape)
         end
      else
         # Skip boundary vertices
      end
   end


   # Plot connections between cell centres and edge midpoints
   # for i=1:nCells
   #    x=findall(x->x>0.5,B̄[i,:])
   #    for j in x
   #       plot!([cellPositions[i,1],edgeMidpoints[j,1]],[cellPositions[i,2],edgeMidpoints[j,2]],color=:red,linewidth=3)
   #    end
   # end

   # Plot connections between neighbouring cell centres
   # xs = zeros(2)
   # ys = zeros(2)
   # for i=1:nEdges
   #    x=findall(x->x>0.5,B̄[:,i])
   #    if size(x)[1] > 1
   #       xs[1] = cellPositions[x[1],1]
   #       xs[2] = cellPositions[x[2],1]
   #       ys[1] = cellPositions[x[1],2]
   #       ys[2] = cellPositions[x[2],2]
   #       plot!(xs,ys,color=:red,alpha=edgeDots[i],linewidth=3)
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

   savefig("output/$folderName/plot$(@sprintf("%03d",outputCount)).png")

end

export visualise

end
