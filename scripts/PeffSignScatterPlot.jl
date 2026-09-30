# Script to make scatter plots of Peff along boundary cells to back claims of sign 

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


# Date string to define the folder we will save to: 
dateString = "26-09-10-10-49-26"
# Vector of paths to data we want to calculate from:
 # List of paths to JLD2 files IN ORDER from (O) to (V)
jld2pathVec = ["data/multipleRuns/26-09-10-10-49-26/(O)/(I)_equilibriumPhase.jld2", #(O)
                "data/multipleRuns/26-09-10-10-49-26/(P)/(NEW 2)_equilibriumPhase.jld2",#(P)
                "data/multipleRuns/26-09-10-10-49-26/(Q)/26-09-12-18-25-13_nCells=1273_Λ_AA=-0.1_Λ_AB=-0.2_Λ_BB=-0.3_β=0.0_γ=0.05/frameData/systemData099.jld2",#(Q)
                "data/multipleRuns/26-09-10-10-49-26/(R)/(III)_equilibriumPhase.jld2",#(R)
                "data/multipleRuns/26-09-10-10-49-26/(S)/26-09-13-03-46-48_nCells=1169_Λ_AA=-0.3_Λ_AB=-0.2_Λ_BB=-0.3_β=0.0_γ=0.05/frameData/systemData099.jld2",#(S)
                "data/multipleRuns/26-09-10-10-49-26/(T)/(V)_equilibriumPhase.jld2",#(T)
                "data/multipleRuns/26-09-10-10-49-26/(U)/(VIII)_equilibriumPhase.jld2",#(U)
                "data/multipleRuns/26-09-10-10-49-26/(V)/(NEW 3)_equilibriumPhase.jld2",]#(V)

# decide which property we would like to scatter 
plotPeffOnBoundary = false
plotξsAcrossMonolayer = true
plotPeffAcrossMonolayer = false

# Orange family — boundary_A relates to exterior
exteriorColor  = RGB(255/255, 178/255, 102/255)   # light orange (existing)
boundaryAColor = RGB(204/255, 102/255,   0/255)   # darker/burnt orange

# Blue family — boundary_B relates to interior
interiorColor  = RGB(102/255, 178/255, 255/255)   # light blue (existing)
boundaryBColor = RGB(  0/255,  76/255, 153/255)   # darker/navy blue


if plotPeffOnBoundary

    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    fig = Figure(size=(1800,600))

    # Initialise figure: 
    grid = fig[1,1] = GridLayout()
    ax = Axis(grid[1,1],
        xticks = (1:8,["(O)","(P)","(Q)","(R)","(S)","(T)","(U)","(V)"]),
        xlabel = "Simulation label",
        ylabel = "ξᵢ")

    parameterLabelVec = ["(O)","(P)","(Q)","(R)","(S)","(T)","(U)","(V)"]

    for (simIdx, jld2pathString) in enumerate(jld2pathVec)

        parameterSetLabel = parameterLabelVec[simIdx]

        R = load(jld2pathString,"R")
        params = load(jld2pathString,"params")
        matrices = load(jld2pathString,"matrices")

        @unpack A, B, C, edgeLabels, cellLabels, P_effs, ξs = matrices
        @unpack nCells = params

        # Initialise vectors for each group: 
        boundaryPeffs_A = []
        boundaryPeffs_B = []

        # Fill in values along the composite boundary: 
        interfaceBoundaryEdges = findall(x -> x==2, edgeLabels)
        interfaceBoundaryCells = Int[]
        for j in interfaceBoundaryEdges
            incidentCells = findall(x->x!=0, @view B[:,j])
            push!(interfaceBoundaryCells, incidentCells...)
        end
        unique!(interfaceBoundaryCells)
        for i in interfaceBoundaryCells
            if cellLabels[i] == 0
                push!(boundaryPeffs_A,P_effs[i])
            elseif cellLabels[i] == 1
                push!(boundaryPeffs_B,P_effs[i])
            end
        end

        # Plot this simulation's three groups immediately, dodged side by side
        boxplot!(ax, fill(simIdx, length(boundaryPeffs_B)), boundaryPeffs_B,
            dodge = fill(2, length(boundaryPeffs_B)), n_dodge = 2, color =interiorColor, width = 0.7)
        boxplot!(ax, fill(simIdx, length(boundaryPeffs_A)), boundaryPeffs_A,
            dodge = fill(3, length(boundaryPeffs_A)), n_dodge = 2, color =exteriorColor, width = 0.7)

    end
    display(fig)
end


if plotξsAcrossMonolayer

    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    fig = Figure(size=(1200,600))

    # Initialise figure: 
    grid = fig[1,1] = GridLayout()
    ax = Axis(grid[1,1],
        xticks = (1:8,["(O)","(P)","(Q)","(R)","(S)","(T)","(U)","(V)"]),
        xlabel = "Simulation label",
        ylabel = "ξᵢ")

    parameterLabelVec = ["(O)","(P)","(Q)","(R)","(S)","(T)","(U)","(V)"]
    # jld2pathVec = jld2pathVec[[1;3:7]] # exclude simulations (P) and (V)

    for (simIdx, jld2pathString) in enumerate(jld2pathVec)

        parameterSetLabel = parameterLabelVec[simIdx]

        R = load(jld2pathString,"R")
        params = load(jld2pathString,"params")
        matrices = load(jld2pathString,"matrices")

        @unpack A, B, C, edgeLabels, cellLabels, P_effs, ξs = matrices
        @unpack nCells = params

        # Initialise vectors for each group: 
        boundaryξs_A = []
        boundaryξs_B = []
        interiorξs = []
        exteriorξs = []

        # Fill in values along the composite boundary: 
        interfaceBoundaryEdges = findall(x -> x==2, edgeLabels)
        interfaceBoundaryCells = Int[]
        for j in interfaceBoundaryEdges
            incidentCells = findall(x->x!=0, @view B[:,j])
            push!(interfaceBoundaryCells, incidentCells...)
        end
        unique!(interfaceBoundaryCells)
        for i in interfaceBoundaryCells
            if cellLabels[i]==0
                push!(boundaryξs_A, ξs[i])
            elseif cellLabels[i]==1
                push!(boundaryξs_B, ξs[i])
            end
            
        end

        # Fill in values on the interior/exterior 
        for i in 1:nCells
            if cellLabels[i] == 0 && !(i in interfaceBoundaryCells)
                push!(exteriorξs, ξs[i])
            elseif cellLabels[i] == 1 && !(i in interfaceBoundaryCells)
                push!(interiorξs, ξs[i])
            end
        end

        
        # Plot this simulation's three groups immediately, dodged side by side
        boxplot!(ax, fill(simIdx, length(boundaryξs_A)), boundaryξs_A,
                dodge = fill(1, length(boundaryξs_A)), n_dodge = 4,
                color = boundaryAColor,
                width = 0.7)
        boxplot!(ax, fill(simIdx, length(exteriorξs)), exteriorξs,
            dodge = fill(2, length(exteriorξs)), n_dodge = 4, color =exteriorColor, width = 0.7)

        boxplot!(ax, fill(simIdx, length(boundaryξs_B)), boundaryξs_B,
            dodge = fill(3, length(boundaryξs_B)), n_dodge = 4, color=boundaryBColor, width = 0.7)
        if !isempty(interiorξs)
            boxplot!(ax, fill(simIdx, length(interiorξs)), interiorξs,
            dodge = fill(4, length(interiorξs)), n_dodge = 4, color =interiorColor, width = 0.7)
        end
        
    end
    save(datadir("multipleRuns",dateString,"ξsGlobalPlot.png"),fig)



end

if plotPeffAcrossMonolayer

    colors = [:blue, :green, :red]
    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    fig = Figure(size=(1200,600))

    # Initialise figure: 
    grid = fig[1,1] = GridLayout()
    ax = Axis(grid[1,1],
        xticks = (1:8,["(O)","(P)","(Q)","(R)","(S)","(T)","(U)","(V)"]),
        xlabel = "Simulation label",
        ylabel = "Peff")

    parameterLabelVec = ["(O)","(P)","(Q)","(R)","(S)","(T)","(U)","(V)"]

    for (simIdx, jld2pathString) in enumerate(jld2pathVec)

        parameterSetLabel = parameterLabelVec[simIdx]

        R = load(jld2pathString,"R")
        params = load(jld2pathString,"params")
        matrices = load(jld2pathString,"matrices")

        @unpack A, B, C, edgeLabels, cellLabels, P_effs, ξs = matrices
        @unpack nCells = params

        # Initialise vectors for each group: 
        boundaryPeffsA = []
        boundaryPeffsB = []
        interiorPeffs = []
        exteriorPeffs = []

        # Fill in values along the composite boundary: 
        interfaceBoundaryEdges = findall(x -> x==2, edgeLabels)
        interfaceBoundaryCells = Int[]
        for j in interfaceBoundaryEdges
            incidentCells = findall(x->x!=0, @view B[:,j])
            push!(interfaceBoundaryCells, incidentCells...)
        end
        unique!(interfaceBoundaryCells)
         for i in interfaceBoundaryCells
            if cellLabels[i]==0
                push!(boundaryPeffsA, P_effs[i])
            elseif cellLabels[i]==1
                push!(boundaryPeffsB, P_effs[i])
            end
            
        end

        # Fill in values on the interior/exterior 
        for i in 1:nCells
            if cellLabels[i] == 0 && !(i in interfaceBoundaryCells)
                push!(exteriorPeffs, P_effs[i])
            elseif cellLabels[i] == 1 && !(i in interfaceBoundaryCells)
                push!(interiorPeffs, P_effs[i])
            end
        end

        # Plot this simulation's three groups immediately, dodged side by side
        boxplot!(ax,fill(simIdx,length(boundaryPeffsA)),boundaryPeffsA,
            dodge = fill(1,length(boundaryPeffsA)),n_dodge = 4,color = boundaryAColor, width=0.7)
        boxplot!(ax, fill(simIdx, length(exteriorPeffs)), exteriorPeffs,
            dodge = fill(2, length(exteriorPeffs)), n_dodge = 4, color =exteriorColor, width = 0.7)
        boxplot!(ax,fill(simIdx,length(boundaryPeffsB)),boundaryPeffsB,
            dodge = fill(3,length(boundaryPeffsB)),n_dodge = 4,color = boundaryBColor, width=0.7)
        if !isempty(interiorPeffs)
            boxplot!(ax, fill(simIdx, length(interiorPeffs)), interiorPeffs,
                dodge = fill(4, length(interiorPeffs)), n_dodge = 4, color =interiorColor, width = 0.7)
        end
        
    end
    save(datadir("multipleRuns",dateString,"PeffsGlobalPlot.png"),fig)
end