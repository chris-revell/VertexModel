using DrWatson
println("Activating environment...")
@quickactivate "VertexModel"
using Revise
using FromFile
println("Loading VertexModel.jl...")
@from "src/VertexModel.jl" using VertexModel
println("Loading test parameters...")
includet("scripts/testParameters.jl")
