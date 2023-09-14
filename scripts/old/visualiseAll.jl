# Import Julia packages
using DrWatson
@quickactivate
using Revise
using LinearAlgebra
using DelimitedFiles
using SparseArrays
using StaticArrays
using CairoMakie
using UnPack
using GeometryBasics
using Random
using Colors
using JLD2
using Printf

# Local modules
# includet("$(projectdir())/src/VertexModelContainers.jl"); using .VertexModelContainers
includet("$(projectdir())/scripts/analysisFunctions/functions.jl")


function visualiseFrame!(dataDirectory,params,matrices,i,t,fig,ax1,ax2,mov,centralCell)
   @unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,outputTotal,realCycleTime,t1Threshold = params
   @unpack A,B,C,R,F,B,Bᵀ,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,externalF,boundaryVertices = matrices
   t+=dt
   empty!(ax1)
   ax1.title="Monolayer in real space"
   # Plot cells
   cellPolygons = makeCellPolygons(params,matrices)
   if plotCells == 1
      for i=1:nCells
         poly!(ax1,cellPolygons[i],color=(getRandomColor(i),0.5))
      end
   end
   scatter!(ax1,Point2f.([cellPositions[centralCell]]),color=:red,markersize=12)
   # annotations!(ax1,"Cell $centralCell", Point2f.([cellPositions[centralCell]]),color=:red)
   # Scatter vertices
   if scatterVertices == 1
      scatter!(ax1,Point2f.(R),color=:green)
      annotations!(ax1,string.(collect(1:length(R))), Point2f.(R),color=:green)
   end
   # Scatter edge midpoints
   if scatterEdges == 1
      scatter!(ax1,Point2f.(edgeMidpoints),color=:blue)
      annotations!(ax1,string.(collect(1:length(edgeMidpoints))), Point2f.(edgeMidpoints),color=:blue)
   end
   # Scatter cell positions
   if scatterCells == 1
      scatter!(ax1,Point2f.(cellPositions),color=:red)
      annotations!(ax1,string.(collect(1:length(cellPositions))), Point2f.(cellPositions),color=:red)
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
            if B[c,i] < 0
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
      arrows!(ax1,xs,us,color=colours,arrowsize=25,linewidth=5)
   end
   if plotForces == 1
      # Plot resultant forces on vertices (excluding external pressure)
      arrows!(ax1,Point2f.(R),Vec2f.(sum(F,dims=2)),color=:green)
      # Plot resultant forces on cells
      arrows!(ax1,Point2f.(cellPositions),Vec2f.(sum(F,dims=1)),color=:red)
   end
   empty!(ax2)
   ax2.title = "Cell $centralCell force space"
   cellNeighbourMatrix = B*Bᵀ
   dropzeros!(cellNeighbourMatrix)
   neighbouringCells = findall(!iszero,cellNeighbourMatrix[centralCell,:])
   # Find and sort all vertices around cells neighbouring centralCell
   cellVerticesDict = Dict()
   for c in neighbouringCells
       # Find vertices around cell
       cellVertices = findall(x->x!=0,C[c,:])
       # Find angles of vertices around cell
       vertexAngles = zeros(size(cellVertices))
       for (k,v) in enumerate(cellVertices)
           vertexAngles[k] = atan((R[v].-cellPositions[c])...)
       end
       # Sort vertices around cell by polar angle
       cellVertices .= cellVertices[sortperm(vertexAngles)]
       # Store sorted cell vertices for this cell
       cellVerticesDict[c] = cellVertices
   end
   centralCellVertices = findall(x->x!=0,C[centralCell,:])
   centralVertexAngles = zeros(size(centralCellVertices))
   for (k,v) in enumerate(centralCellVertices)
      centralVertexAngles[k] = atan((R[v].-cellPositions[centralCell])...)
   end
   m = minimum(centralVertexAngles)
   centralVertexAngles .-= m
   centralCellVertices .= centralCellVertices[sortperm(centralVertexAngles)]
   # Sort cells neighbouring centralCell by angle
   setdiff!(neighbouringCells,[centralCell]) # Remove centralCell from neighbours list
   neighbourAngles = zeros(length(neighbouringCells))
   for (i,c) in enumerate(neighbouringCells)
       neighbourAngles[i] = atan((cellPositions[c].-cellPositions[centralCell])...)
   end
   neighbourAngles .+= (2π-m)
   neighbourAngles = neighbourAngles.%(2π)
   neighbouringCells .= neighbouringCells[sortperm(neighbourAngles)]
   # Draw force network
   startPosition = @SVector [0.0,0.0]
   for (i,v) in enumerate(cellVerticesDict[centralCell])
       arrows!(ax2,Point2f.([startPosition]),Vec2f.([ϵ*F[v,centralCell]]),linewidth=6,arrowsize=24,color=(getRandomColor(centralCell),0.75))
       annotateForceSpace == 1 ? annotations!(ax2,string.([v]),Point2f.([startPosition.+ϵ*F[v,centralCell]./2.0]),color=(getRandomColor(centralCell),0.75)) : nothing
       startPosition = startPosition + ϵ*F[v,centralCell]
       H = Array{SVector{2,Float64}}(undef,length(cellVerticesDict[neighbouringCells[i]])+1)
       cellForces = SVector{2, Float64}[]
       # Circular permutation of vertices to ensure vertex v is the first index
       # in the ordered cellVertices list around cell neighbouringCells[i]
       index = findall(x->x==v, cellVerticesDict[neighbouringCells[i]])
       cellVertices = circshift(cellVerticesDict[neighbouringCells[i]],1-index[1])
       H[1] = startPosition
       for (j,cv) in enumerate(cellVertices)
           push!(cellForces,+ϵ*F[cv,neighbouringCells[i]])
           H[j+1] = H[j]+cellForces[end]
       end
       annotateForceSpace == 1 ? annotations!(ax2,string.(cellVertices),(Point2f.(H[1:end-1])+Vec2f.(cellForces)./2.0),color=(getRandomColor(neighbouringCells[i]),0.75)) : nothing
       arrows!(ax2,Point2f.(H),Vec2f.(cellForces),color=(getRandomColor(neighbouringCells[i]),0.75),linewidth=6,arrowsize=24)
   end
   recordframe!(mov)
   save("$dataDirectory/frames/frame$(@sprintf("%03d", i)).png",fig)
end


# #dataDirectory = "data/sims/2022-02-28-19-30-22"
dataDirs = ["data/figure7/$x" for x in readdir("data/figure7/") if isdir("data/figure7/$x")]

plotCells         = 1
plotEdges         = 0
scatterEdges      = 0
scatterVertices   = 0
scatterCells      = 0
plotForces        = 0
annotateForceSpace= 0

centralCell = 1

dataDirectory = dataDirs[3]

# for dataDirectory in dataDirs
   fig = Figure(resolution=(2000,3000),fontsize = 32)
   fig = Figure(resolution=(2000,1000))
   grid = fig[1,1] = GridLayout()
   ax1 = Axis(grid[1,1],aspect=DataAspect())
   ax2 = Axis(grid[1,2],aspect=DataAspect())
   Label(fig[2,1,Bottom()],"Movie showing the monolayer from Figure 7(g-i) in real space alongside \nthe network of rotated forces acting at all vertices of cell $centralCell.\n Cell $centralCell highlighted in the monolayer with a red dot.",fontsize = 32)


   resize_to_layout!(fig)
   hidedecorations!(ax1)
   hidespines!(ax1)
   hidedecorations!(ax2)
   hidespines!(ax2)

   # colsize!(grid, 1, Aspect(1, 1))
   # colsize!(grid, 2, Aspect(1, 1))
   # colsize!(fig.layout, 3, Aspect(1, 1.5))
   # colsize!(fig.layout, 1, Aspect(1, 1))
   # rowsize!(fig.layout, 1, Aspect(1, 0.5))
   rowsize!(fig.layout, 2, Aspect(1, 0.02))
   # colgap!(fig.layout,Relative(0.0))
   # rowgap!(fig.layout,Relative(0.01))

   mov = VideoStream(fig, framerate=5)
   t=0.0


   for i=0:100
      conditionsDict    = load("$dataDirectory/frames/params$(@sprintf("%03d", i)).jld2")
      matricesDict = load("$dataDirectory/frames/matrices$(@sprintf("%03d", i)).jld2")
      visualiseFrame!(dataDirectory,conditionsDict["params"],matricesDict["matrices"],i,t,fig,ax1,ax2,mov,centralCell)
   end

   save("$dataDirectory/Movie1.mp4",mov)
# end