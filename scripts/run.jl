# Julia packages
using DrWatson
@quickactivate
using Revise
include("$(projectdir())/scripts/testParameters.jl")
using VertexModel
vertexModel("seven",4.0*realTimetMax,realCycleTime,γ,λ,viscousTimeScale,dt,preferredArea,pressureExternal,outputTotal,t1Threshold,outputToggle)
