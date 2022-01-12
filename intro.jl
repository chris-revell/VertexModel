using DrWatson
@quickactivate "VertexModel"
using Revise
if "src/" in LOAD_PATH ? nothing : push!(LOAD_PATH,"src/")
using VertexModel
includet("scripts/testParameters.jl")
println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())
"""
)
