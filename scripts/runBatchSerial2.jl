using DrWatson
@quickactivate
using VertexModel

realTimetMax     = 1.5*86400.0
γ                = [0.15,0.2,0.1]
L₀               = [3.0,0.75,0.5]
outputTotal      = 100
outputToggle     = 1
frameDataToggle  = 0
frameImageToggle = 0
printToggle      = 0
videoToggle      = 1
subFolder        = "examples"

for i=1:3
    vertexModel(realTimetMax = realTimetMax, γ = γ[i], L₀ = L₀[i], outputTotal = outputTotal, outputToggle = outputToggle, frameDataToggle = frameDataToggle, frameImageToggle = frameImageToggle, printToggle = printToggle, videoToggle = videoToggle, subFolder = subFolder)
    vertexModel(realTimetMax = realTimetMax, γ = γ[i], L₀ = L₀[i], outputTotal = outputTotal, outputToggle = outputToggle, frameDataToggle = frameDataToggle, frameImageToggle = frameImageToggle, printToggle = printToggle, videoToggle = videoToggle, subFolder = subFolder)
end
