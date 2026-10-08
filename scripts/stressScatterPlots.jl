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
# Vector of paths to data we want to calculate from:
 # List of paths to JLD2 files IN ORDER from (O) to (V)
jld2pathVec = ["data/multipleRuns/26-09-10-10-49-26/(O)/26-09-12-14-42-03_nCells=1273_Λ_AA=-0.2_Λ_AB=-0.2_Λ_BB=-0.2_β=0.0_γ=0.05/frameData/systemData099.jld2", #(O)
                "data/multipleRuns/26-09-10-10-49-26/(P)/26-09-12-13-53-32_nCells=1273_Λ_AA=-0.1_Λ_AB=-0.2_Λ_BB=-0.25_β=0.0_γ=0.05/frameData/systemData099.jld2",#(P)
                "data/multipleRuns/26-09-10-10-49-26/(Q)/26-09-12-18-25-13_nCells=1273_Λ_AA=-0.1_Λ_AB=-0.2_Λ_BB=-0.3_β=0.0_γ=0.05/frameData/systemData099.jld2",#(Q)
                "data/multipleRuns/26-09-10-10-49-26/(R)/26-09-12-13-39-48_nCells=1273_Λ_AA=-0.15_Λ_AB=-0.2_Λ_BB=-0.3_β=0.0_γ=0.05/frameData/systemData099.jld2",#(R)
                "data/multipleRuns/26-09-10-10-49-26/(S)/26-09-13-03-46-48_nCells=1169_Λ_AA=-0.3_Λ_AB=-0.2_Λ_BB=-0.3_β=0.0_γ=0.05/frameData/systemData099.jld2",#(S)
                "data/multipleRuns/26-09-10-10-49-26/(T)/26-09-12-14-02-50_nCells=1273_Λ_AA=-0.3_Λ_AB=-0.2_Λ_BB=-0.15_β=0.0_γ=0.05/frameData/systemData099.jld2",#(T)
                "data/multipleRuns/26-09-10-10-49-26/(U)/26-09-12-14-08-47_nCells=1273_Λ_AA=-0.3_Λ_AB=-0.2_Λ_BB=-0.1_β=0.0_γ=0.05/frameData/systemData099.jld2",#(U)
                "data/multipleRuns/26-09-10-10-49-26/(V)/26-09-12-13-57-33_nCells=1273_Λ_AA=-0.25_Λ_AB=-0.2_Λ_BB=-0.1_β=0.0_γ=0.05/frameData/systemData099.jld2",]#(V)


parameterLabelVec = ["(O)","(P)","(Q)","(R)","(S)","(T)","(U)","(V)"]

plotξs = true

for (simInd, jld2pathString) in enumerate(jld2pathVec)

    

    parameterSetLabel = parameterLabelVec[simInd]

    isdir(datadir("multipdleRuns",dateString,"CoupeStressScatterPlots",parameterSetLabel)) ? nothing : mkpath(datadir("multipleRuns",dateString,"CoupeStressScatterPlots",parameterSetLabel)) 

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
            ξs,
            edgeLengths = matrices

    # First compute couple stresses: 
    h = hNetwork(R,A,B,F)
    coupleStresses = -curlᵛ(R,A,B,h)

    interfaceBoundaryEdges = findall(x -> x==2,edgeLabels)

    # Initialise vectors to store edge and pressure differences along interfaceBoundaryEdges
    PeffDiffVec = zeros(Float64,length(interfaceBoundaryEdges))
    PeffDiffVecTimesLength = zeros(Float64,length(interfaceBoundaryEdges))
    CoupleStressDiffVec = zeros(Float64,length(interfaceBoundaryEdges))
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

        for cell in cells 
            # Always take first entry as the B-cell (interior on comp boundary)
            if cellLabels[cell] == 1
                incidentCells[1]=cell
            elseif cellLabels[cell] == 0
                incidentCells[2]=cell
            end
        end
        orderAroundBCell , ~ = OrderAroundCell.orderAroundCell(matrices,incidentCells[1])

        trailingVertices = findall(x->x!=0, @view(A[j_oldInd,:]))
        # Check which order these vertices appear in going clockwise around cell B:
        positionVert1 = findfirst(x -> x == trailingVertices[1], orderAroundBCell)
        positionVert2 = findfirst(x -> x == trailingVertices[2], orderAroundBCell)

        # Need to account for the fact this is a circular array - check whether the index is at the start/end
        n = length(orderAroundBCell)
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
        PeffDiffVecTimesLength[j_newInd] =PeffDiffVec[j_newInd] * edgeLengths[j_oldInd]
        CoupleStressDiffVec[j_newInd] = coupleStresses[incidentVerts[2]] - coupleStresses[incidentVerts[1]]
        ξDiffVec[j_newInd] = (ξs[incidentCells[2]] - ξs[incidentCells[1]])

        

    end

    # Initialise scatter plot
    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    fig = Figure(size=(600,600))

    # Initialise a figure for tracking sum of P_effsA_i: 
    grid = fig[1,1] = GridLayout()
    ax = Axis(grid[1,1],aspect=1)
    ax.title = "$parameterSetLabel Effective pressure difference against couple stress difference across interface edges"
    ax.xlabel = "ΔP_eff*lⱼ"
    ax.ylabel = "Δ{CURLh}ₖ"

    scatter!(ax, PeffDiffVecTimesLength, CoupleStressDiffVec, color=:blue, markersize=5)

    X = [ones(length(PeffDiffVec)) PeffDiffVecTimesLength]
    β = X \ CoupleStressDiffVec   # least squares solutions
    intercept, slope = β
    xs = range(extrema(PeffDiffVecTimesLength)..., length=200)
    lines!(ax, xs, intercept .+ slope .* xs, color=:black, linewidth=2)

    save(datadir("multipleRuns",dateString,"CoupeStressScatterPlots",parameterSetLabel, "PeffCSDiffScatterPlot.png"), fig)

    if plotξs

        # Initialise scatter plot
        set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
        fig = Figure(size=(600,600))

        # Initialise a figure for tracking sum of P_effsA_i: 
        grid = fig[1,1] = GridLayout()
        ax = Axis(grid[1,1],aspect=1)
        ax.title = "$parameterSetLabel Sheer stress difference against couple stress difference across interface edges"
        ax.xlabel = "Δξ"
        ax.ylabel = "Δ{CURLh}ₖ"

        scatter!(ax, ξDiffVec, CoupleStressDiffVec, color=:blue, markersize=5)

        X = [ones(length(ξDiffVec)) ξDiffVec]
        β = X \ CoupleStressDiffVec   # least squares solution
        intercept, slope = β
        xs = range(extrema(ξDiffVec)..., length=200)
        lines!(ax, xs, intercept .+ slope .* xs, color=:black, linewidth=2)

        save(datadir("multipleRuns",dateString,"CoupeStressScatterPlots",parameterSetLabel, "ξCSDiffScatterPlot.png"),fig)


    end

    

    

end
