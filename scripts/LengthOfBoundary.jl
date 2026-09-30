# Script to run over edges of a system and calculate the length of the interfacial boundary 

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


# Pick the equilibrated system: 
dateString = "26-09-10-10-49-26"
parameterSetLabel = "(O)" 

!isdir(datadir("multipleRuns", dateString)) ? mkpath(datadir("multipleRuns", dateString)) : nothing 

# dataDict = load(datadir("multipleRuns",dateString, parameterSetLabel, "$(parameterSetLabel)_equilibriumPhase.jld2");
#                     typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer" => MatricesContainer, 
#                                 "VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer" => ParametersContainer
#                     )
#                 )

# dataDict = load("C:/Users/z28439ct/Documents/GitHub/VertexModel/data/multipleRuns/26-09-10-10-49-26/(I)/26-09-12-14-42-03_nCells=1273_Λ_AA=-0.2_Λ_AB=-0.2_Λ_BB=-0.2_β=0.0_γ=0.05/frameData/systemData099.jld2";
#                 typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer" => MatricesContainer, 
#                             "VertexModel...VertexModelContainers.jl.VertexModelContainers.ParametersContainer" => ParametersContainer))

jld2pathVec = ["data/multipleRuns/26-09-10-10-49-26/(O)/(I)_equilibriumPhase.jld2", #(O)
                "data/multipleRuns/26-09-10-10-49-26/(P)/(NEW 2)_equilibriumPhase.jld2",#(P)
                "data/multipleRuns/26-09-10-10-49-26/(Q)/26-09-12-18-25-13_nCells=1273_Λ_AA=-0.1_Λ_AB=-0.2_Λ_BB=-0.3_β=0.0_γ=0.05/frameData/systemData099.jld2",#(Q)
                "data/multipleRuns/26-09-10-10-49-26/(R)/(III)_equilibriumPhase.jld2",#(R)
                "data/multipleRuns/26-09-10-10-49-26/(S)/26-09-13-03-46-48_nCells=1169_Λ_AA=-0.3_Λ_AB=-0.2_Λ_BB=-0.3_β=0.0_γ=0.05/frameData/systemData099.jld2",#(S)
                "data/multipleRuns/26-09-10-10-49-26/(T)/(V)_equilibriumPhase.jld2",#(T)
                "data/multipleRuns/26-09-10-10-49-26/(U)/(VIII)_equilibriumPhase.jld2",#(U)
                "data/multipleRuns/26-09-10-10-49-26/(V)/(NEW 3)_equilibriumPhase.jld2",#(V)
                "data/multipleRuns/26-09-10-10-49-26/Symmetric S/systemData048.jld2",#Symmetric S
                "data/sims/charlie-free-boundaries/26-09-23-15-45-29_nCells=213_Λ_AA=-0.35_Λ_AB=-0.2_Λ_BB=-0.3_β=0.0_γ=0.05/frameData/systemData015.jld2"]# Extreme case

boundaryLengthVec = []
boundaryLengthDivArea = []
surfaceTensionVec = []

for (simInd, jld2pathString) in enumerate(jld2pathVec)

    isdir(datadir("multipleRuns",dateString,"CoupeStressScatterPlots",parameterSetLabel)) ? nothing : mkpath(datadir("multipleRuns",dateString,"CoupeStressScatterPlots",parameterSetLabel)) 

    R = load(jld2pathString,"R")
    params = load(jld2pathString,"params")
    matrices = load(jld2pathString,"matrices")

    @unpack edgeLabels,
        edgeLengths,
        cellAreas,
        cellLabels = matrices
    @unpack Λ_AA,
        Λ_AB,
        Λ_BB,
        nCells = params
    

     interfaceBoundaryEdges = findall(x -> x==2,edgeLabels)
     interfaceLength = 0
     for j in interfaceBoundaryEdges
        interfaceLength += edgeLengths[j]
     end

     BCellAreaSum = 0
     for i = 1:nCells
        if cellLabels[i] == 1
            BCellAreaSum += cellAreas[i]
        end
     end
     surfaceTension = Λ_AB - (Λ_AA + Λ_BB)/2



     push!(surfaceTensionVec,surfaceTension)
     push!(boundaryLengthVec,interfaceLength)
     push!(boundaryLengthDivArea,interfaceLength/sqrt(BCellAreaSum))


end

println(surfaceTensionVec)
# Initialise scatter plot
set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
fig = Figure(size=(1200,600))

# Initialise a figure for tracking sum of P_effsA_i: 
grid = fig[1,1] = GridLayout()
ax = Axis(grid[1,1],aspect=1)
ax2 = Axis(grid[1,2],aspect=1)
ax.title = "Interface length against heterotypic interfacial tension"
ax.xlabel = "γ_(AB) "
ax.ylabel = "(Σⱼlⱼ), j on heterotypic interface"
ax2.xlabel = "γ_(AB) "
ax2.ylabel = "(Σⱼlⱼ)/sqrt(ΣᵢAᵢ), j on heterotypic interface i in B"

scatter!(ax,surfaceTensionVec,boundaryLengthVec,color=:blue, markersize=5)
scatter!(ax2,surfaceTensionVec,boundaryLengthDivArea,color=:blue,markersize=5)


# X = [ones(length(surfaceTensionVec)) surfaceTensionVec]
# # β = X \ boundaryLengthVec   # least squares solution
# intercept, slope = β
# xs = range(extrema(surfaceTensionVec)..., length=200)
# lines!(ax, xs, intercept .+ slope .* xs, color=:black, linewidth=2)

save(datadir("multipleRuns",dateString,"InterfaceLengthAgainstSurfaceTension.png"),fig)