using DrWatson
using Revise
using FromFile

println("Activating environment...")
@quickactivate "VertexModel"

println("Loading VertexModel.jl...")
@from "src/VertexModel.jl" using VertexModel

println("Loading test parameters...")
includet("scripts/testParameters.jl")

println("Warming up function...")
vertexModel(initialSystem,10.0,realCycleTime,γ[1],L₀[1],A₀,viscousTimeScale,0.1,pressureExternal,peripheralTension,t1Threshold,1,0,0;subFolder="")