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
parameterSetLabel = "(III)"
# Path to data we want to calculate from:
jld2pathString = "/Users/user/The University of Manchester Dropbox/Charlotte Taylor Barca/JULIA/VertexModel/data/sims/charlie-free-boundaries/Symmetric simulations/NEW 2/26-09-02-18-08-20_nCells=469_Λ_AA=-0.1_Λ_AB=-0.2_Λ_BB=-0.25_β=0.0_γ=0.05/frameData/systemData036.jld2"

R = load(jld2pathString,"R")
params = load(jld2pathString,"params")
matrices = load(jld2pathString,"matrices")

@unpack A,
        B,
        C,
        F,
        edgeLabels,
        cellLabels,
        P_effs = matrices

# First compute couple stresses: 
h = hNetwork(R,A,B,F)
coupleStresses = -curlᵛ(R,A,B,h)

interfaceBoundaryEdges = findall(x -> x==2,edgeLabels)

# Initialise vectors to store edge and pressure differences along interfaceBoundaryEdges
PeffDiffVec = zeros(Float64,length(interfaceBoundaryEdges))
CoupleStressDiffVec = zeros(Float64,length(interfaceBoundaryEdges))

for edge in enumerate(interfaceBoundaryEdges)
    # store new and old indices: 
    j_newInd = edge[1]
    j_oldInd = edge[2]

    incidentCells = zeros(Int64,2)
    incidentVerts = zeros(Int64,2)

    # Find incident cells and then store them in order A-B 
    cells = findall(x -> x!=0, @view B[:,j_oldInd])
    orderAroundCell_B = CircularArray{Int64}
    for cell in cells
        if cellLabels[cell] == 0 # cell A
            incidentCells[1] = cell
            # orderAroundCell_A = OrderAroundCell.orderAroundCell(matrices,cell)
        else
            incidentCells[2] = cell
            # Find a list of vertices going clockwise around cell B, we want to take from left to right when looking out from cell B.
            orderAroundCell_B, ~ = OrderAroundCell.orderAroundCell(matrices,cell)
        end
    end

    trailingVertices = findall(x->x!=0, @view(A[j_oldInd,:]))
    # Check which order these vertices appear in going clockwise around cell B:
    positionVert1 = findall(x -> x == trailingVertices[1], orderAroundCell_B)
    positionVert2 = findall(x -> x == trailingVertices[2], orderAroundCell_B)

    if positionVert1 < positionVert2
        incidentVerts[1] = trailingVertices[1]
        incidentVerts[2] = trailingVertices[2]
    else
        incidentVerts[1] = trailingVertices[2]
        incidentVerts[2] = trailingVertices[1]
    end

    PeffDiffVec[j_newInd] = P_effs[incidentCells[1]] - P_effs[incidentCells[2]]
    CoupleStressDiffVec[j_newInd] = coupleStresses[incidentVerts[1]] - coupleStresses[incidentVerts[2]]
    
end

# Initialise scatter plot
set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
fig = Figure(size=(600,600))

# Initialise a figure for tracking sum of P_effsA_i: 
grid = fig[1,1] = GridLayout()
ax = Axis(grid[1,1],aspect=1)
ax.title = "Effective pressure difference against couple stress differnce across interface edges"
ax.xlabel = "ΔP_eff"
ax.ylabel = "Δ{CURLh}ₖ"

scatter!(ax, PeffDiffVec, CoupleStressDiffVec, color=:blue, markersize=5)
display(fig)

# save(datadir("multipleRuns","dateString","parameterSetLabel", "stressDiffScatterPlot.png"), fig)
