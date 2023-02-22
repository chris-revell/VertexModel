using DrWatson
using FromFile
@quickactivate "VertexModel"
@from "$(projectdir())/src/VertexModel.jl" using VertexModel
includet("$(projectdir())/scripts/testParameters.jl")
# vertexModel(initialSystem,realTimetMax*5.0,realCycleTime,γ[1],L₀[1],A₀,viscousTimeScale,dt,pressureExternal,peripheralTension*0.0,t1Threshold*5,outputTotal,outputToggle,plotToggle;subFolder="annealing")
vertexModel(initialSystem,realTimetMax,realCycleTime,γ[2],L₀[2],A₀,viscousTimeScale,pressureExternal,peripheralTension*0.0,t1Threshold,outputTotal,outputToggle,plotToggle;subFolder="test")
# vertexModel(initialSystem,realTimetMax*5.0,realCycleTime,γ[3],L₀[3],A₀,viscousTimeScale,dt,pressureExternal,peripheralTension*0.0,t1Threshold*5,outputTotal,outputToggle,plotToggle;subFolder="annealing")