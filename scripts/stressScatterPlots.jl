# Script to make scatter plots of the change in effective pressure against the couple stress difference on edges 
# along the population interface from equilibrium. 

# We take care in the manner in which we compute differences (e.g. from A-to-B, from left-to-right)

using DrWatson
using DiscreteCalculus
using CairoMakie
using StaticArrays
using VertexModel
using LinearAlgebra
using SparseArrays
using Random
using Colors 
using JLD2
using Dates
using FromFile
using InvertedIndices
using LaTeXStrings
using Dates 
using CircularArrays
using OrdinaryDiffEq
using Printf
using DiffEqCallbacks

include(srcdir("OrderAroundCell.jl")); using .OrderAroundCell

# Date string to define the folder we will save to: 
dateString = "26-09-10-10-49-26"
parameterSetLabel = "(XIII)"
# Path to data we want to calculate from:
jld2pathString = "data/multipleRuns/26-09-10-10-49-26/(XIII)/26-09-13-03-46-48_nCells=1169_Λ_AA=-0.3_Λ_AB=-0.2_Λ_BB=-0.3_β=0.0_γ=0.05/frameData/systemData099.jld2"


R = load(jld2pathString,"R")
params = load(jld2pathString,"params")
matrices = load(jld2pathString,"matrices")

@unpack A,
        B,
        C,
        F,
        edgeLabels,
        cellLabels,
        P_effs,
        ξs = matrices

# First compute couple stresses: 
h = hNetwork(R,A,B,F)
coupleStresses = -curlᵛ(R,A,B,h)

interfaceBoundaryEdges = findall(x -> x==2,edgeLabels)

# Initialise vectors to store edge and pressure differences along interfaceBoundaryEdges
PeffDiffVec = zeros(Float64,length(interfaceBoundaryEdges))
CoupleStressDiffVec = zeros(Float64,length(interfaceBoundaryEdges))
CoupleStressDiffVec2 = zeros(Float64,length(interfaceBoundaryEdges))
ξDiffVec = zeros(Float64,length(interfaceBoundaryEdges))

for edge in enumerate(interfaceBoundaryEdges)
    # store new and old indices: 
    j_newInd = edge[1]
    j_oldInd = edge[2]

    incidentCells = zeros(Int64,2)
    incidentVerts = zeros(Int64,2)

    incidentCells2 = zeros(Int64,2)
    incidentVerts2 = zeros(Int64,2)

    # Find incident cells and then store them in order A-B 
    cells = findall(x -> x!=0, @view B[:,j_oldInd])
    orderAroundLowerPeffCell = CircularArray{Int64}
    orderAroundLowerPeffCell2 = CircularArray{Int64}

    # Check which of the cells has smaller Peff foor ordering: 
    if P_effs[cells[1]] < P_effs[cells[2]]
        incidentCells[1] = cells[1]
        incidentCells[2] = cells[2]
        orderAroundLowerPeffCell, ~ = OrderAroundCell.orderAroundCell(matrices,cells[1])
    else
        incidentCells[1] = cells[2]
        incidentCells[2] = cells[1]
        orderAroundLowerPeffCell, ~ = OrderAroundCell.orderAroundCell(matrices,cells[2])
    end

    if ξs[cells[1]] < ξs[cells[2]]
        incidentCells2[1] = cells[1]
        incidentCells2[2] = cells[2]
        orderAroundLowerPeffCell2, ~ = OrderAroundCell.orderAroundCell(matrices,cells[1])
    else
        incidentCells2[1] = cells[2]
        incidentCells2[2] = cells[1]
        orderAroundLowerPeffCell2, ~ = OrderAroundCell.orderAroundCell(matrices,cells[2])
    end

    trailingVertices = findall(x->x!=0, @view(A[j_oldInd,:]))
    # Check which order these vertices appear in going clockwise around cell B:
    positionVert1 = findfirst(x -> x == trailingVertices[1], orderAroundLowerPeffCell)
    positionVert2 = findfirst(x -> x == trailingVertices[2], orderAroundLowerPeffCell)

    # Need to account for the fact this is a circular array - check whether the index is at the start/end
    n = length(orderAroundLowerPeffCell)
    if mod(positionVert2 - positionVert1, n) == 1
        # forward (CW) traversal goes trailingVertices[1] -> trailingVertices[2]
        incidentVerts[1] = trailingVertices[1]
        incidentVerts[2] = trailingVertices[2]
    elseif mod(positionVert1 - positionVert2, n) == 1
         # forward (CW) traversal goes trailingVertices[2] -> trailingVertices[1]
        incidentVerts[1] = trailingVertices[2]
        incidentVerts[2] = trailingVertices[1]
    else
        error("Vertices $(trailingVertices) are not adjacent in orderAroundLowerPeffCell — check edge/cell correspondence.")
    end
    

    PeffDiffVec[j_newInd] = P_effs[incidentCells[2]] - P_effs[incidentCells[1]]
    ξDiffVec[j_newInd] = ξs[incidentCells[2]] - ξs[incidentCells[1]]
    CoupleStressDiffVec[j_newInd] = coupleStresses[incidentVerts[2]] - coupleStresses[incidentVerts[1]]
    
    positionVert1b = findfirst(x -> x == trailingVertices[1], orderAroundLowerPeffCell2)
    positionVert2b = findfirst(x -> x == trailingVertices[2], orderAroundLowerPeffCell2)
    n2 = length(orderAroundLowerPeffCell2)

    incidentVerts2 = zeros(Int64,2)
    if mod(positionVert2b - positionVert1b, n2) == 1
        incidentVerts2[1] = trailingVertices[1]
        incidentVerts2[2] = trailingVertices[2]
    elseif mod(positionVert1b - positionVert2b, n2) == 1
        incidentVerts2[1] = trailingVertices[2]
        incidentVerts2[2] = trailingVertices[1]
    else
        error("Vertices $(trailingVertices) are not adjacent in orderAroundLowerPeffCell2.")
    end

    CoupleStressDiffVec2[j_newInd] = coupleStresses[incidentVerts2[2]] - coupleStresses[incidentVerts2[1]]


end

# Initialise scatter plot
set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
fig = Figure(size=(1200,600))

# Initialise a figure for tracking sum of P_effsA_i: 
grid = fig[1,1] = GridLayout()
ax = Axis(grid[1,1],aspect=1)
ax2 = Axis(grid[1,2],aspect=1)
ax.title = "$parameterSetLabel Effective pressure difference against couple stress difference across interface edges"
ax2.title = "$parameterSetLabel Shear stress difference against couple stress difference across interface edges"
ax.xlabel = "ΔP_eff"
ax2.xlabel = "Δξ"
ax.ylabel = "Δ{CURLh}ₖ"

scatter!(ax, PeffDiffVec, CoupleStressDiffVec, color=:blue, markersize=5)
scatter!(ax2, ξDiffVec, CoupleStressDiffVec2, color=:blue, markersize=5)

save(datadir("multipleRuns",dateString,parameterSetLabel, "stressDiffScatterPlot.png"), fig)
