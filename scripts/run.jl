using DrWatson
using FromFile
@quickactivate "VertexModel"
@from "$(projectdir())/src/VertexModel.jl" using VertexModel
includet("$(projectdir())/scripts/testParameters.jl")
vertexModel()