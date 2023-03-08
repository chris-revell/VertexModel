using DrWatson
using Revise
using FromFile

println("Activating environment...")
@quickactivate "VertexModel"

println("Loading VertexModel.jl...")
@from "src/VertexModel.jl" using VertexModel

println("Loading test parameters...")
includet("scripts/testParameters.jl")