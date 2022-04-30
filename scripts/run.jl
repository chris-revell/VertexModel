# Julia packages
using DrWatson
@quickactivate
using Revise
include("$(projectdir())/scripts/testParameters.jl")
using VertexModel

vertexModel(initialSystem,realTimetMax,realCycleTime,γ,L₀,A₀,viscousTimeScale,dt,pressureExternal,t1Threshold,outputTotal,outputToggle,subFolder="")
