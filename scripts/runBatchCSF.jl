# Julia packages
using DrWatson
@quickactivate
# using Revise
# using Base.Threads
using VertexModel

initialSystem    = "seven"#"data/sims/AlexPaperParameters/2022-04-22-15-12-34/"
realTimetMax     = 4.0*86400.0 # = 24 hours
realCycleTime    = 86400.0 # = 24 hours
γ                = 0.2#0.15#0.1#
L₀               = 0.25#3.0#0.5#         0.27#0.03#-0.9#-0.1#-0.3
viscousTimeScale = 20.0
dt               = 0.5
A₀               = 1.0
pressureExternal = 0.0
outputTotal      = 100
t1Threshold      = 0.01
outputToggle     = 1
plotToggle       = 0

@threads for seed=1:8
    vertexModel(initialSystem,realTimetMax,realCycleTime,γ,L₀,A₀,viscousTimeScale,dt,pressureExternal,t1Threshold,outputTotal,outputToggle,0;subFolder="Test")
end
