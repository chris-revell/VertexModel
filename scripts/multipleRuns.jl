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

deNovoSystem = true 

if deNovoSystem

    integ0 = vertexModel(initialSystem = "new",
    boundaryType = "free",
    cellLayout = "random",
    nRows = 3,
    nCycles = 10,
    realCycleTime = 86400.0, # From Megan's data, using the division rate 0.15/min 
    realTimetMax = nCycles*realCycleTime,
    viscousTimeScale = 1000.0,
    β = 0.1,
    divisionToggle = 0,
    nBlasThreads = 1,
    subFolder = "",
    outputTotal = 100,
    outputToggle = 1,
    frameDataToggle = 1,
    frameImageToggle = 1,
    printToggle = 1,
    videoToggle = 1,
    plotCells = 1,
    scatterEdges = 0,
    scatterVertices = 0,
    scatterCells = 0,
    plotXis = 1,
    plotStresses = 1,
    plotForces = 0,
    plotEdgeMidpointLinks = 0,
    randomSeed = 0,
    energyModel = "quadratic2pops",
    vertexWeighting = 1,
    noiseWeighting = 1,
    R_in = spzeros(2),
    A_in = spzeros(2),
    B_in = spzeros(2), 
    L_x = 20,
    L_y = 20,
    Λ_AA = -0.2, 
    Λ_AB = -0.2,
    Λ_BB = -0.2,
    Λ_AE = -0.2, 
    Λ_BE = -0.2,
    Area_A_ratio = 0.5,
    t1timeGap = 1e-0,
    spiky = true,
    desiredNumCells = 100,
    plotOrientations = 1,
    edgeToAblate = [],
    clusterWidth = 5, # the radius of the central B cluster in cell number, in the case initialSystem = "symmetric" 
    termSteadyState = false, # flag to determine whether simulation terminates once it reaches steady state 
    randomDivision = true, # flag to determine whether division process is random or not (i.e., cell cycle times are uniform)
)

end
