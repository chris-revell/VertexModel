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
jld2pathVec = ["data/multipleRuns/26-09-10-10-49-26/(III)/(III)_equilibriumPhase.jld2"]


# decide which property we would like to scatter 
plotPeffOnBoundary = true 

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


    for jld2pathString in jld2pathVec
        # Unpack params and matrices for this path: 
        R = load(jld2pathString,"R")
        params = load(jld2pathString,"params")
        matrices = load(jld2pathString,"matrices")

        @unpack A,
                B,
                C,
                edgeLabels,
                cellLabels,
                P_effs,
                ξs = matrices

        @unpack Λ_AA,
                Λ_BB = params

        interfaceBoundaryEdges = findall(x -> x==2,edgeLabels)
        interfacePeffVecA = []
        interfacePeffVecB = []
        for j in interfaceBoundaryEdges
            incidentCells = findall(x->x!=0,@view B[:,j])
            for cell in incidentCells
                if cellLabels[cell] == 0
                    push!(interfacePeffVecA,P_effs[cell])
                else
                    push!(interfacePeffVecB,P_effs[cell])
                end
            end
        end
        unique!(interfacePeffVecA)
        unique!(interfacePeffVecB)


        Λ_AA_vec = Λ_AA*ones(length(interfacePeffVecA))
        Λ_BB_vec = Λ_BB*ones(length(interfacePeffVecB))
        
        scatter!(Λ_AAPeffBoundaryAx,Λ_AA_vec,interfacePeffVecA,color=:blue,markersize=5)
        scatter!(Λ_BBPeffBoundaryAx,Λ_BB_vec,interfacePeffVecB,color=:red,markersize=5)
    end
    
    display(fig)
end