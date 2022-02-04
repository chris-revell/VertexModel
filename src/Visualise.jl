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
using StaticArrays
using SparseArrays

function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end

@views function visualise(t,fig,ax,ax2,mov,params,matrices)

   plotCells         = 1
   plotEdges         = 0
   scatterEdges      = 0
   scatterVertices   = 0
   scatterCells      = 1
   plotForces        = 0

   @unpack R,A,B,C,cellPositions,edgeTangents,edgeMidpoints,F,ϵ = matrices
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
   # For each edge, use A incidence matrix to find corresponding vertices x, and plot line between x[1] and x[2]
   if plotEdges == 1
      xs = Point2f[]
      us = Vec2f[]
      colours = Tuple{Symbol, Float64}[]
      for c=1:nCells
         es = findall(x->x!=0,B[c,:])
         for i in es
            vs = findall(x->x!=0,A[i,:])
            colour=:black
            # Use B to set colour of edge depending on whether it runs with or against the orientation of the cell face
            if matrices.B[c,i] < 0
               colour = (:red,0.25)
            else
               colour = (:blue,0.25)
            end
            # Use A to set direction of arrow along edge
            if A[i,vs[1]] < 0
               push!(xs, Point2f(R[vs[1]]))
               push!(us, Vec2f(edgeTangents[i]))
               push!(colours, colour)
            else
               push!(xs, Point2f(R[vs[2]]))
               push!(us, Vec2f(edgeTangents[i]))
               push!(colours, colour)
            end
         end
      end
      arrows!(ax,xs,us,color=colours,arrowsize=25,linewidth=5)
   end

   if plotForces == 1
      # Plot resultant forces on vertices (excluding external pressure)
      arrows!(ax,Point2f.(R),Vec2f.(sum(F,dims=2)),color=:green)
      # Plot resultant forces on cells
      arrows!(ax,Point2f.(matrices.cellPositions),Vec2f.(sum(F,dims=1)),color=:red)
   end


   empty!(ax2)

   ax2.title = "Cell 36 force space"

   centralCell = 36

   cellNeighbourMatrix = matrices.B*matrices.Bᵀ
   dropzeros!(cellNeighbourMatrix)
   neighbouringCells = findall(!iszero,cellNeighbourMatrix[centralCell,:])

   # Find and sort all vertices around cells neighbouring centralCell
   cellVerticesDict = Dict()
   for c in neighbouringCells
       # Find vertices around cell
       cellVertices = findall(x->x!=0,matrices.C[c,:])
       # Find angles of vertices around cell
       vertexAngles = zeros(size(cellVertices))
       for (k,v) in enumerate(cellVertices)
           vertexAngles[k] = atan((R[v].-matrices.cellPositions[c])...)
       end
       # Sort vertices around cell by polar angle
       cellVertices .= cellVertices[sortperm(vertexAngles)]
       # Store sorted cell vertices for this cell
       cellVerticesDict[c] = cellVertices
   end

   # Sort cells neighbouring centralCell by angle
   setdiff!(neighbouringCells,[centralCell]) # Remove centralCell from neighbours list
   neighbourAngles = zeros(length(neighbouringCells))
   for (i,c) in enumerate(neighbouringCells)
       neighbourAngles[i] = atan((cellPositions[c].-cellPositions[centralCell])...)
   end
   neighbouringCells .= neighbouringCells[sortperm(neighbourAngles)]

   # Draw force network
   startPosition = @SVector [0.0,0.0]
   for (i,v) in enumerate(cellVerticesDict[centralCell])
       arrows!(ax2,Point2f.([startPosition]),Vec2f.([ϵ*F[v,centralCell]]),linewidth=4,arrowsize=16,color=(getRandomColor(centralCell),0.75))
       annotations!(ax2,string.([v]),Point2f.([startPosition.+ϵ*F[v,centralCell]./2.0]),color=(getRandomColor(centralCell),0.75))
       startPosition = startPosition + ϵ*F[v,centralCell]
       H = Array{SVector{2,Float64}}(undef,length(cellVerticesDict[neighbouringCells[i]])+1)
       cellForces = SVector{2, Float64}[]
       # Circular permutation of vertices to ensure vertex v is the first index
       # in the ordered cellVertices list around cell neighbouringCells[i]
       index = findall(x->x==v, cellVerticesDict[neighbouringCells[i]])
       cellVertices = circshift(cellVerticesDict[neighbouringCells[i]],1-index[1])
       H[1] = startPosition
       for (j,cv) in enumerate(cellVertices)
           push!(cellForces,ϵ*F[cv,neighbouringCells[i]])
           H[j+1] = H[j]+cellForces[end]
       end
       annotations!(ax2,string.(cellVertices),(Point2f.(H[1:end-1])+Vec2f.(cellForces)./2.0),color=(getRandomColor(neighbouringCells[i]),0.75))
       arrows!(ax2,Point2f.(H),Vec2f.(cellForces),color=(getRandomColor(neighbouringCells[i]),0.75),linewidth=4,arrowsize=16)
   end

   recordframe!(mov)

   return nothing

end

export visualise

end
