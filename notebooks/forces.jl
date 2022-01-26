
using DrWatson
@quickactivate
using Revise
display(projectdir())
using LinearAlgebra
using DelimitedFiles
using SparseArrays
using StaticArrays
using CairoMakie
using UnPack
using GeometryBasics
using DelimitedFiles
using Random
using Colors
using JLD2

# Local modules
includet("$(projectdir())/src/TopologyChange.jl"); using .TopologyChange
includet("$(projectdir())/src/CreateRunDirectory.jl"); using .CreateRunDirectory
includet("$(projectdir())/src/Visualise.jl"); using .Visualise
includet("$(projectdir())/src/Initialise.jl"); using .Initialise
includet("$(projectdir())/src/Iterate.jl"); using .Iterate
includet("$(projectdir())/src/SpatialData.jl"); using .SpatialData

ϵ                 = @SMatrix [  # Rotation matrix setting orientation of cell faces
    0.0 -1.0
    1.0 0.0
]

function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end

# initialSystem = "data/sims/2022-01-24-20-29-36"
initialSystem = "/Users/christopher/Dropbox (The University of Manchester)/Other/VertexModelStuff/2022-01-25-12-48-54"
#initialSystem      = "$(projectdir())/data/sims/2022-01-20-14-52-51"

conditionsDict    = load("$initialSystem/params.jld2")

@unpack γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea = conditionsDict["params"]
pressureExternal   = 0.1
outputTotal        = 100
realCycleTime      = 24.0*60*60
t1Threshold        = 0.1

params,matrices = initialise(initialSystem,realTimetMax,γ,λ,preferredArea,pressureExternal,dt,viscousTimeScale,outputTotal,t1Threshold,realCycleTime)

topologyChange!(matrices)

spatialData!(matrices.R,params,matrices)

matrixF = Array{SVector{2,Float64}}(undef,params.nVerts,params.nCells)

function calculateForce!(params,matrices,matrixF)
    @unpack R,A,B,Ā,B̄,cellTensions,cellPressures,edgeLengths,edgeTangents,F,ϵ,boundaryVertices = matrices
    @unpack nVerts,nCells,nEdges,pressureExternal = params
    fill!(matrixF,@SVector zeros(2))
    for k=1:nVerts
        for i=1:nCells
            for j=1:nEdges
                matrixF[k,i] -= 0.5*cellPressures[i]*B[i,j]*Ā[j,k]*(ϵ*edgeTangents[j]) - cellTensions[i]*B̄[i,j]*A[j,k]*edgeTangents[j]/edgeLengths[j]
            end
        end
    end
end

calculateForce!(params,matrices,matrixF)

#%% This section plots the full system as fig1

fig1 = Figure(resolution=(1000,1000))
grid1 = fig1[1,1] = GridLayout()
ax1 = Axis(grid1[1,1],aspect=DataAspect(),)
hidedecorations!(ax1)
hidespines!(ax1)
for i=1:params.nCells
  cellVertices = findall(x->x!=0,matrices.C[i,:])
  vertexAngles = zeros(size(cellVertices))
  for (k,v) in enumerate(cellVertices)
     vertexAngles[k] = atan((matrices.R[v].-matrices.cellPositions[i])...)
  end
  cellVertices .= cellVertices[sortperm(vertexAngles)]
  poly!(ax1,Point2f.(matrices.R[cellVertices]),color=(getRandomColor(i),0.25))
end
scatter!(ax1,Point2f.(matrices.R),alpha=0.5,color=:black,size=0.5)
scatter!(ax1,Point2f.(matrices.cellPositions),color=:red)
annotations!(ax1,string.(collect(1:params.nCells)),Point2f.(matrices.cellPositions),color=:red)
display(fig1)

reinterpretedR = reinterpret(reshape, Float64, matrices.R)


# Display force vectors on each vertex from neighbouring cells
# for v=1:params.nVerts
#     #reinterpretedFv = reinterpret(reshape, Float64, matrixF[v,:])
#     nonZeroIndices = findall(!iszero,norm.(matrixF[v,:]))
#     cellForce = @SVector zeros(2)
#     for (n,f) in enumerate(nonZeroIndices)
#         vertexForce += matrixF[v,f]
#     end
#     display(vertexForce)
#     # arrows!(ax1,[reinterpretedR[1,v]],[reinterpretedR[2,v]],[3.0*vertexForce[1]],[3.0*vertexForce[2]],linecolor=:blue,alpha=0.3)
# end

# Display force vectors on each cell from surrounding vertices
for c=1:params.nCells
    nonZeroIndices = findall(!iszero,norm.(matrixF[:,c]))
    cellForce = @SVector zeros(2)
    for (n,f) in enumerate(nonZeroIndices)
        cellForce += matrixF[f,c]
    end
    #display(vertexForce)
    arrows!(ax1,[matrices.cellPositions[c][1]],[matrices.cellPositions[c][2]],[cellForce[1]],[cellForce[2]],linecolor=:red)
end

#%% This section plots, for a given cell, the forces from each neighbouring cell on each vertex, and the force on the cell from each vertex

# Work with cell 31
cell=31

fig2 = Figure(resolution=(1000,1000))
grid2 = fig2[1,1] = GridLayout()
ax2 = Axis(grid2[1,1],aspect=DataAspect(),)
hidedecorations!(ax2)
hidespines!(ax2)

# Find all vertices around the cell
cellVertices = findall(x->x!=0,matrices.C[cell,:])
# Find all cells around each vertex of original cell
vertexCells = Int64[]
for v in cellVertices
    append!(vertexCells,findall(x->x!=0,matrices.C[:,v]))
end
# Remove duplicates to obtain list of neighbouring cells

neighbouringCells = unique(vertexCells)
for c in neighbouringCells
    if c != cell
        cellVertices = findall(x->x!=0,matrices.C[c,:])
        vertexAngles = zeros(size(cellVertices))
        for (k,v) in enumerate(cellVertices)
            vertexAngles[k] = atan((matrices.R[v].-matrices.cellPositions[c])...)
        end
        cellVertices .= cellVertices[sortperm(vertexAngles)]
        poly!(ax2,Point2f.(matrices.R[cellVertices]),color=(getRandomColor(c),0.25))
    end
end

cellVertices = findall(x->x!=0,matrices.C[cell,:])
vertexAngles = zeros(size(cellVertices))
for (k,v) in enumerate(cellVertices)
    vertexAngles[k] = atan((matrices.R[v].-matrices.cellPositions[cell])...)
end
cellVertices .= cellVertices[sortperm(vertexAngles)]
poly!(ax2,Point2f.(matrices.R[cellVertices]),color=(getRandomColor(cell),0.25))

scatter!(ax2,Point2f.(matrices.R[cellVertices]),color=:black)
annotations!(ax2,string.(cellVertices),Point2f.(matrices.R[cellVertices]),color=:black)

scatter!(ax2,Point2f.(matrices.cellPositions[neighbouringCells]),color=:red)
annotations!(ax2,string.(neighbouringCells),Point2f.(matrices.cellPositions[neighbouringCells]),color=:red)

# Display force vectors on each vertex from neighbouring cells
for (i,v) in enumerate(cellVertices)
    nonZeroIndices = findall(!iszero,norm.(matrixF[v,:]))
    vertexForce = @SVector zeros(2)
    for f in nonZeroIndices
        rotatedForce = ϵ*matrixF[v,f]
        vertexForce += rotatedForce
        arrows!(ax2,Point2f.([matrices.R[v]]),Vec2f.([rotatedForce]),color=:blue)
        annotations!(ax2,["c$f"],Point2f.([(matrices.R[v].+rotatedForce)]),color=:blue)
    end
    arrows!(ax2,Point2f.([matrices.R[v]]),Vec2f.([vertexForce]),color=:red)
end

# Display force vectors on each cell from surrounding vertices
nonZeroIndices = findall(!iszero,norm.(matrixF[:,cell]))
cellForceResultant = @SVector zeros(2)
for f in nonZeroIndices
    rotatedForce = ϵ*matrixF[f,cell]
    cellForceResultant += rotatedForce
    arrows!(ax2,Point2f.([matrices.cellPositions[cell]]),Vec2f.([rotatedForce]),color=:blue)
    annotations!(ax2,["v$f"],Point2f.([(matrices.cellPositions[cell].+rotatedForce)]),color=:blue)
end
arrows!(ax2,Point2f.([matrices.cellPositions[cell]]),Vec2f.([cellForce]),linecolor=:red)

display(fig2)

#%%
