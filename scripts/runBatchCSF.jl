# Julia packages
using DrWatson
@quickactivate
using Revise
using Base.Threads
include("$(projectdir())/scripts/testParameters.jl")
using VertexModel

for seed=1:6
    vertexModel(initialSystem,realTimetMax,realCycleTime,γ,L₀,A₀,viscousTimeScale,dt,pressureExternal,t1Threshold,outputTotal,outputToggle,0;subFolder="Test")
end
