# Julia packages
using DrWatson
@quickactivate
using Revise
using Base.Threads
include("$(projectdir())/scripts/testParameters.jl")
using VertexModel

@threads for seed=1:6
    vertexModel(seed,initialSystem,realTimetMax,realCycleTime,γ,λ,viscousTimeScale,dt,A₀,pressureExternal,outputTotal,t1Threshold,outputToggle)
end
