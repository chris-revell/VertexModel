using DrWatson
@quickactivate
using VertexModel
using DifferentialEquations

nRuns = 1

initialSystem = "large"
realTimetMax = 2.0 * 86400.0
realCycleTime = 86400.0
γ = [0.15, 0.2, 0.1]
L₀ = [3.0, 0.75, 0.5]
A₀ = 1.0
viscousTimeScale = 20.0
pressureExternal = 0.0
peripheralTension = 0.0
t1Threshold = 0.05
solver = Tsit5()
nBlasThreads = 1
subFolder = ""
outputTotal = 100
outputToggle = 1
frameDataToggle = 1
frameImageToggle = 1
printToggle = 1
videoToggle = 1
plotCells = 1
scatterEdges = 0
scatterVertices = 0
scatterCells = 0
plotForces = 0
setRandomSeed = 0

for _ = 1:nRuns
    vertexModel(
        initialSystem = initialSystem,
        realTimetMax = realTimetMax,
        realCycleTime = realCycleTime,
        γ = γ[1],
        L₀ = L₀[1],
        A₀ = A₀,
        viscousTimeScale = viscousTimeScale,
        pressureExternal = pressureExternal,
        peripheralTension = peripheralTension,
        t1Threshold = t1Threshold,
        solver = solver,
        nBlasThreads = nBlasThreads,
        subFolder = subFolder,
        outputTotal = outputTotal,
        outputToggle = outputToggle,
        frameDataToggle = frameDataToggle,
        frameImageToggle = frameImageToggle,
        printToggle = printToggle,
        videoToggle = videoToggle,
        plotCells = plotCells,
        scatterEdges = scatterEdges,
        scatterVertices = scatterVertices,
        scatterCells = scatterCells,
        plotForces = plotForces,
        setRandomSeed = setRandomSeed
    )
end
