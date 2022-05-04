# Julia packages
using DrWatson
@quickactivate

initialSystem    = "seven"
realTimetMax     = 5.0*86400.0
realCycleTime    = 86400.0
γ                = [0.2,0.15,0.1]
L₀               = [0.75,3.0,0.5]
viscousTimeScale = 20.0
dt               = 0.5
A₀               = 1.0
pressureExternal = 0.1
outputTotal      = 100
t1Threshold      = 0.01
outputToggle     = 1

using VertexModel

for i=1:3
    vertexModel(initialSystem,realTimetMax,realCycleTime,γ[1],L₀[1],A₀,viscousTimeScale,dt,pressureExternal,t1Threshold,outputTotal,outputToggle,subFolder="Figure7")
end
