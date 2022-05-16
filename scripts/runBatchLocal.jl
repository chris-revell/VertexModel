# Julia packages
using DrWatson
@quickactivate

initialSystem    = "seven"
realTimetMax     = 4.0*86400.0
realCycleTime    = 86400.0
γ                = 0.20
L₀               = 0.75
viscousTimeScale = 20.0
dt               = 0.5
A₀               = 1.0
pressureExternal = 0.2
outputTotal      = 100
t1Threshold      = 0.1
outputToggle     = 1
plotToggle       = 1
folderLabel      = "InitialParameterRuns"

using VertexModel

for seed=1:10
    vertexModel(initialSystem,realTimetMax,realCycleTime,γ,L₀,A₀,viscousTimeScale,dt,pressureExternal,t1Threshold,outputTotal,outputToggle,plotToggle;subFolder=folderLabel)
end
