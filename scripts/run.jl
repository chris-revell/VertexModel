# Julia packages
using DrWatson
@quickactivate
using Revise
include("$(projectdir())/scripts/testParameters.jl")
using VertexModel
vertexModel(initialSystem,realTimetMax,realCycleTime,γ,λ,viscousTimeScale,dt,preferredArea,pressureExternal,outputTotal,t1Threshold,outputToggle)
