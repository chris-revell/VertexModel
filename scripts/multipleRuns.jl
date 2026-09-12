# Script to run multiple instances of vertex model for different purposes 

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

deNovoSystem = false # starting a new monolayer from scratch
fixedSeed = 12 # Fixed across parameter conditions 

runningFromStage = 1 # stage 1: running the 

if deNovoSystem

    # This is what we will run to get the starting system for all parameter cases

    integ0 = vertexModel(initialSystem ="32-cell",
                    nRows = 3,
                    nCycles = 6,
                    realCycleTime = 86400.0, 
                    viscousTimeScale = 1000.0,
                    β = 0.1,
                    divisionToggle = 1,
                    outputToggle = 1,
                    frameDataToggle = 0,
                    frameImageToggle = 1,
                    videoToggle = 1,
                    energyModel = "quadratic2pops",
                    Λ_AA = -0.2, 
                    Λ_AB = -0.2,
                    Λ_BB = -0.2,
                    Λ_AE = -0.2, 
                    Λ_BE = -0.2,
                    t1timeGap = 1e-0,
                    spiky = true,
                    termSteadyState = false, # flag to determine whether simulation terminates once it reaches steady state 
                    randomDivision = true, # flag to determine whether division process is random or not (i.e., cell cycle times are uniform)
                    randomSeed = fixedSeed,
                )

    R = reinterpret(SVector{2,Float64}, integ0.u)
    (params, matrices) = integ0.p

    dateString = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))"
    !isdir(datadir("multipleRuns", dateString)) ? mkpath(datadir("multipleRuns", dateString)) : nothing 
    jldsave(datadir("multipleRuns", dateString, "$(dateString)_InitialSystem.jld2"); 
                R,
                params,
                matrices
            )

else

    dateString = "26-09-10-10-49-26"
    dataDict = load(datadir("multipleRuns", dateString, "$(dateString)_InitialSystem.jld2");
                    typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer" => MatricesContainer, 
                                "VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer" => ParametersContainer
                    )
                )

    # Import system data
    @unpack R, params, matrices = dataDict

    

end


# Set desired values for line tensions of this run: 
params.Λ_AA = -0.1
params.Λ_BB = -0.1
params.Λ_AE = -0.1
params.Λ_BE = -0.1

parameterSetLabel = "(NEW)"
!isdir(datadir("multipleRuns", dateString, parameterSetLabel)) ? mkpath(datadir("multipleRuns", dateString,parameterSetLabel)) : nothing 

γ = params.γ

integ1 = vertexModel(initialSystem = "argument",
                    nCycles = 3.3,
                    realCycleTime = 86400.0, 
                    viscousTimeScale = 1000.0,
                    β = 0.1,
                    divisionToggle = 1,
                    outputToggle = 1,
                    frameDataToggle = 1,
                    frameImageToggle = 1,
                    printToggle = 1,
                    videoToggle = 1,
                    γ=γ,
                    R_in = R,
                    A_in = matrices.A,
                    B_in = matrices.B,
                    Λ_AA = params.Λ_AA, 
                    Λ_AB = params.Λ_AB,
                    Λ_BB = params.Λ_BB,
                    Λ_AE = params.Λ_AE, 
                    Λ_BE = params.Λ_BE,
                    termSteadyState = false, # flag to determine whether simulation terminates once it reaches steady state 
                    randomDivision = false, # flag to determine whether division process is random or not (i.e., cell cycle times are uniform)
                    randomSeed = fixedSeed,
                )

R = reinterpret(SVector{2,Float64}, integ1.u)
(params, matrices) = integ1.p

jldsave(datadir("multipleRuns", dateString, parameterSetLabel, "$(parameterSetLabel)_growthPhase.jld2"); matrices,params,R)

integ2 = vertexModel(initialSystem = datadir("multipleRuns", dateString, parameterSetLabel, "$(parameterSetLabel)_growthPhase.jld2"),
                    realCycleTime = 86400.0, 
                    viscousTimeScale = 1000.0,
                    nCycles = 2,
                    β = 0.1,
                    divisionToggle = 0,
                    frameDataToggle = 1,
                    frameImageToggle = 1,
                    printToggle = 1,
                    videoToggle = 1,
                    γ=γ,
                    R_in = R,
                    A_in = matrices.A,
                    B_in = matrices.B,
                    Λ_AA = params.Λ_AA, 
                    Λ_AB = params.Λ_AB,
                    Λ_BB = params.Λ_BB,
                    Λ_AE = params.Λ_AE, 
                    Λ_BE = params.Λ_BE,
                    termSteadyState = false, # flag to determine whether simulation terminates once it reaches steady state 
                    randomSeed = fixedSeed,
                )

R = reinterpret(SVector{2,Float64}, integ2.u)
(params, matrices) = integ2.p

jldsave(datadir("multipleRuns", dateString, parameterSetLabel, "$(parameterSetLabel)_noGrowthPhase.jld2"); matrices,params,R)

integ3 = vertexModel(initialSystem = datadir("multipleRuns", dateString, parameterSetLabel, "$(parameterSetLabel)_noGrowthPhase.jld2"),
                    realCycleTime = 86400.0, 
                    viscousTimeScale = 1000.0,
                    β = 0.0,
                    divisionToggle = 0,
                    frameDataToggle = 0,
                    frameImageToggle = 1,
                    printToggle = 1,
                    videoToggle = 0,
                    termSteadyState = true, # flag to determine whether simulation terminates once it reaches steady state 
                    randomSeed = fixedSeed,
                )

jldsave(datadir("multipleRuns", dateString, parameterSetLabel, "$(parameterSetLabel)_equilibriumPhase.jld2"); matrices,params,R)

