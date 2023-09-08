using DrWatson
using FromFile
# @quickactivate "VertexModel"
using VertexModel

include("$(scriptsdir("testParameters.jl"))")
vertexModel(γ=γ[1], L₀=L₀[1], subFolder="examples")