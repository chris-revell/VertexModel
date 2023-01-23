using DrWatson
using FromFile
@quickactivate "VertexModel"
@from "$(projectdir())/src/VertexModel.jl" using VertexModel
includet("$(projectdir())/scripts/testParameters.jl")
vertexModel(initialSystem,10.0,realCycleTime,γ,L₀,A₀,viscousTimeScale,0.1,pressureExternal,peripheralTension,t1Threshold,1,0,0;subFolder="")
vertexModel(initialSystem,realTimetMax*1.5,realCycleTime,γ,L₀,A₀,viscousTimeScale,dt,pressureExternal,peripheralTension,t1Threshold,outputTotal,outputToggle,plotToggle;subFolder="runs")