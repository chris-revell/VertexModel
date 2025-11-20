# Julia packages
using DrWatson
using FromFile
using OrdinaryDiffEq
using LinearAlgebra
using JLD2
using SparseArrays
using StaticArrays
using CairoMakie
using Printf

# Local modules
@from "CreateRunDirectory.jl" using CreateRunDirectory
@from "Visualise.jl" using Visualise
@from "Initialise.jl" using Initialise
@from "SpatialData.jl" using SpatialData
@from "PlotSetup.jl" using PlotSetup
@from "Model.jl" using Model
@from "T1Transitions.jl" using T1Transitions
@from "TopologyChange.jl" using TopologyChange
@from "Division.jl" using Division
@from "SenseCheck.jl" using SenseCheck
@from "Callbacks.jl" using Callbacks


integ = vertexModel(
    nRows = 9,
    nCycles = 1,
    divisionToggle = 1,
    outputTotal = 1,
    outputToggle = 0,
    frameDataToggle = 0,
    frameImageToggle = 0,
    printToggle = 1,
    plotCells = 0,
    setRandomSeed = 0,
    energyModel = "quadratic",
)



for i = 1:nCells
            poly!(ax, cellPolygons[i], color=(getRandomColor(i), 0.5), strokecolor=(:black, 1.0), strokewidth=2)
        end