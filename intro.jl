push!(LOAD_PATH,"./")
push!(LOAD_PATH,"src")
using DrWatson
println("Activating environment...")
@quickactivate "VertexModel"
using Revise
println("Loading VertexModel.jl...")
#includet("src/VertexModel.jl")
#using .VertexModel
using VertexModel
println("Loading test parameters...")
includet("scripts/testParameters.jl")
