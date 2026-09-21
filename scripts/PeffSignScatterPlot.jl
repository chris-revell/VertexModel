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
parameterSetLabel = "(I)"
# Vector of paths to data we want to calculate from:
jld2pathVec = [] # List of paths to JLD2 files IN ORDER from (O) to (V)

# decide which property we would like to scatter 
plotPeffOnBoundary = false
plotξsAcrossMonolayer = true

if plotPeffOnBoundary

    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    fig = Figure(size=(1200,600))

    # Initialise figure: 
    grid = fig[1,1] = GridLayout()
    Λ_AAPeffBoundaryAx = Axis(grid[1,1],aspect=1)
    Λ_AAPeffBoundaryAx.title = "Values of Peff along the composite boundary against Λ_AA"
    Λ_AAPeffBoundaryAx.xlabel = "Λ_AA"
    Λ_AAPeffBoundaryAx.ylabel = "P_eff on composite boundary"

    Λ_BBPeffBoundaryAx = Axis(grid[1,2],aspect=1)
    Λ_BBPeffBoundaryAx.title = "Values of Peff along the composite boundary against Λ_BB"
    Λ_BBPeffBoundaryAx.xlabel = "Λ_BB"
    Λ_BBPeffBoundaryAx.ylabel = "P_eff on composite boundary"


    Λ_AA_all = Float64[]
    Peff_A_all = Float64[]
    Λ_BB_all = Float64[]
    Peff_B_all = Float64[]

    for jld2pathString in jld2pathVec
        R = load(jld2pathString,"R")
        params = load(jld2pathString,"params")
        matrices = load(jld2pathString,"matrices")

        @unpack A, B, C, edgeLabels, cellLabels, P_effs, ξs = matrices
        @unpack Λ_AA, Λ_BB = params

        interfaceBoundaryEdges = findall(x -> x==2, edgeLabels)
        interfacePeffVecA = Float64[]
        interfacePeffVecB = Float64[]
        for j in interfaceBoundaryEdges
            incidentCells = findall(x->x!=0, @view B[:,j])
            for cell in incidentCells
                if cellLabels[cell] == 0
                    push!(interfacePeffVecA, P_effs[cell])
                else
                    push!(interfacePeffVecB, P_effs[cell])
                end
            end
        end
        unique!(interfacePeffVecA)
        unique!(interfacePeffVecB)

        append!(Λ_AA_all, fill(Λ_AA, length(interfacePeffVecA)))
        append!(Peff_A_all, interfacePeffVecA)
        append!(Λ_BB_all, fill(Λ_BB, length(interfacePeffVecB)))
        append!(Peff_B_all, interfacePeffVecB)
    
    
    end

    boxplot!(Λ_AAPeffBoundaryAx, Λ_AA_all, Peff_A_all,
        width=0.03, show_outliers=true, color=(:blue,0.4), strokecolor=:blue, strokewidth=1)

    boxplot!(Λ_BBPeffBoundaryAx, Λ_BB_all, Peff_B_all,
        width=0.03, show_outliers=true, color=(:red,0.4), strokecolor=:red, strokewidth=1)

    # Overlay median trend line
    using Statistics
    for (Λvec, Peffvec, ax, col) in ((Λ_AA_all, Peff_A_all, Λ_AAPeffBoundaryAx, :blue),
                                    (Λ_BB_all, Peff_B_all, Λ_BBPeffBoundaryAx, :red))
        uniqueΛ = sort(unique(Λvec))
        medians = [median(Peffvec[Λvec .== λ]) for λ in uniqueΛ]
        lines!(ax, uniqueΛ, medians, color=col, linewidth=2)
        scatter!(ax, uniqueΛ, medians, color=col, markersize=8, marker=:diamond)
    end
    
    display(fig)
end

if plotξsAcrossMonolayer

    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    fig = Figure(size=(1200,600))

    # Initialise figure: 
    grid = fig[1,1] = GridLayout()
    ax = Axis(grid[1,1],
        xticks = (1:6,["(O)","(Q)","(R)","(S)","(T)","(U)"]),
        xlabel = "Simulation label",
        ylabel = "ξᵢ")

    parameterLabelVec = ["(O)","(Q)","(R)","(S)","(T)","(U)"]
    jld2pathVec = jld2pathVec[[1;3:7]] # exclude simulations (P) and (V)

    for (simIdx, jld2pathString) in enumerate(jld2pathVec)

        parameterSetLabel = parameterLabelVec[simIdx]

        R = load(jld2pathString,"R")
        params = load(jld2pathString,"params")
        matrices = load(jld2pathString,"matrices")

        @unpack A, B, C, edgeLabels, cellLabels, P_effs, ξs = matrices
        @unpack nCells = params

        # Initialise vectors for each group: 
        boundaryξs = []
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
            push!(boundaryξs, ξs[i])
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
        boxplot!(ax, fill(simIdx, length(boundaryξs)), boundaryξs,
            dodge = fill(1, length(boundaryξs)), color = colors[1], width = 0.7)
        boxplot!(ax, fill(simIdx, length(interiorξs)), interiorξs,
            dodge = fill(2, length(interiorξs)), color = colors[2], width = 0.7)
        boxplot!(ax, fill(simIdx, length(exteriorξs)), exteriorξs,
            dodge = fill(3, length(exteriorξs)), color = colors[3], width = 0.7)

    end



end