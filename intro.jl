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
#
# println(
# """
# Currently active project is: $(projectname())
#
# Path of active project: $(projectdir())
# """
# )
