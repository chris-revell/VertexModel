# Julia packages
using DrWatson
@quickactivate
using Revise
include("$(projectdir())/scripts/testParameters.jl")
using VertexModel
vertexModel("data/sims/2022-02-28-19-30-22",0.5*realTimetMax,realCycleTime,γ,λ,viscousTimeScale,dt,preferredArea,pressureExternal,outputTotal,t1Threshold,outputToggle)
