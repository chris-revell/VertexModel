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
vertexModel(initialSystem,realTimetMax,10.0,γ,L₀,A₀,viscousTimeScale,0.1,pressureExternal,t1Threshold,1,0,0;subFolder="")