
using DrWatson
@quickactivate
using VertexModel

initialSystem    = "seven"
realTimetMax     = 4.0*86400.0
realCycleTime    = 86400.0
γ                = [0.2,0.15,0.1]
L₀               = [0.75,3.0,1.0]
viscousTimeScale = 20.0
dt               = 0.5
A₀               = 1.0
pressureExternal = 0.2
outputTotal      = 100
t1Threshold      = 0.1
outputToggle     = 1
plotToggle       = 1
subFolderName    = "BasicParameterRuns"

for seed=1:6
    vertexModel(initialSystem,realTimetMax,realCycleTime,γ[1],L₀[1],A₀,viscousTimeScale,dt,pressureExternal,t1Threshold,outputTotal,outputToggle,plotToggle;subFolder=subFolderName)
end
# for seed=1:6
#     vertexModel(initialSystem,realTimetMax,realCycleTime,γ[2],L₀[2],A₀,viscousTimeScale,dt,pressureExternal,t1Threshold,outputTotal,outputToggle,plotToggle;subFolder="BatchSerialLocal")
# end
# for seed=1:6
#     vertexModel(initialSystem,realTimetMax,realCycleTime,γ[3],L₀[3],A₀,viscousTimeScale,dt,pressureExternal,t1Threshold,outputTotal,outputToggle,plotToggle;subFolder="BatchSerialLocal")
# end