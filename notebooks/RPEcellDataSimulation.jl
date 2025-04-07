
using LinearAlgebra
using CairoMakie
using StatsBase

@from "$(srcdir("VertexModel.jl"))" using VertexModel
@from "$(srcdir("RotationMatrix.jl"))" using RotationMatrix

folderName = "curved-surface/25-04-04-14-39-13_L₀=0.75_nCells=792_realTimetMax=86400.0_γ=0.2"

R, matrices, params = loadData(datadir("sims",folderName), outputNumber=99)

CairoMakie.activate!()

nBins = 10

cellXY = [norm(a[1:2]) for a in matrices.cellPositions]
cellRadii = params.surfaceRadius.*asin.(cellXY./params.surfaceRadius)
binsEdges = collect(range(minimum(cellRadii), maximum(cellRadii), nBins+1))
binCentres = (binsEdges[1:end-1].+binsEdges[2:end])./2.0

fig = CairoMakie.Figure(size=(1000,1000))

ax0 = Axis(fig[0,1])
# for i=2:length(binsEdges)
#     indices = findall(x->x>=binsEdges[i-1] && x<binsEdges[i], cellRadii)
#     push!(cellAreaMeans, mean(matrices.cellAreas[indices]))
#     push!(cellAreaErrors, sem(matrices.cellAreas[indices]))
# end
scatter!(ax0, cellRadii, matrices.cellL₀s, color=(:blue, 0.1))
# lines!(ax1, binCentres, cellAreaMeans, color=(:black, 1.0))
# errorbars!(ax1, binCentres, cellAreaMeans, cellAreaErrors, color=(:black, 1.0))
ax0.xlabel = "Cell radius from system centre"
ax0.ylabel = "Cell L₀s"

ax1 = Axis(fig[1,1])
cellAreaMeans = Float64[]
cellAreaErrors = Float64[]
for i=2:length(binsEdges)
    indices = findall(x->x>=binsEdges[i-1] && x<binsEdges[i], cellRadii)
    push!(cellAreaMeans, mean(matrices.cellAreas[indices]))
    push!(cellAreaErrors, sem(matrices.cellAreas[indices]))
end
scatter!(ax1, cellRadii, matrices.cellAreas, color=(:blue, 0.1))
lines!(ax1, binCentres, cellAreaMeans, color=(:black, 1.0))
errorbars!(ax1, binCentres, cellAreaMeans, cellAreaErrors, color=(:black, 1.0))
ax1.xlabel = "Cell radius from system centre"
ax1.ylabel = "Cell area"

ax2 = Axis(fig[2,1])
neighbourMatrix = matrices.B*matrices.Bᵀ
cellNeighbourCount = [neighbourMatrix[i,i] for i in 1:params.nCells]
cellNeighbourMeans = Float64[]
cellNeighbourErrors = Float64[]
for i=2:length(binsEdges)
    indices = findall(x->x>=binsEdges[i-1] && x<binsEdges[i], cellRadii)
    push!(cellNeighbourMeans, mean(cellNeighbourCount[indices]))
    push!(cellNeighbourErrors, sem(cellNeighbourCount[indices]))
end
scatter!(ax2, cellRadii, cellNeighbourCount, color=(:blue, 0.1))
lines!(ax2, binCentres, cellNeighbourMeans, color=(:black, 1.0))
errorbars!(ax2, binCentres, cellNeighbourMeans, cellNeighbourErrors, color=(:black, 1.0))
ax2.xlabel = "Cell radius from system centre"
ax2.ylabel = "Neighbour count"

ax3 = Axis(fig[3,1])
cellAspectRatios = []
for i=1:params.nCells
    spokes = [R[kk].-matrices.cellPositions[i] for kk in matrices.cellVertexOrders[i][0:end]]
    cellPerpAxis = matrices.cellPositions[i].-params.surfaceCentre
    crossVec = cellPerpAxis×[1,0,0]
    ϵCoordinates = ϵ(v=crossVec, θ=asin(norm(crossVec)/(norm(cellPerpAxis))))
    rotatedSpokes = [(ϵCoordinates*s)[2:end] for s in spokes]
    cellShapeTensor = sum(rotatedSpokes[2:end].*transpose.(rotatedSpokes[2:end]))./matrices.cellEdgeCount[i]
    # Long and short axis from eigenvectors of shapetensor
    # Put some sort of tolerance that if eigenvalues are approx equal we randomly choose a division orientation, eg circ >0.95
    eigenVals, eigenVecs = eigen(cellShapeTensor) 
    push!(cellAspectRatios, eigenVals[2]/eigenVals[1])
end
cellAspectRatioMeans = Float64[]
cellAspectRatioErrors = Float64[]
for i=2:length(binsEdges)
    indices = findall(x->x>=binsEdges[i-1] && x<binsEdges[i], cellRadii)
    push!(cellAspectRatioMeans, mean(cellAspectRatios[indices]))
    push!(cellAspectRatioErrors, sem(cellAspectRatios[indices]))
end
scatter!(ax3, cellRadii, cellAspectRatios, color=(:blue, 0.1))
lines!(ax3, binCentres, cellAspectRatioMeans, color=(:black, 1.0))
errorbars!(ax3, binCentres, cellAspectRatioMeans, cellAspectRatioErrors, color=(:black, 1.0))
ax3.xlabel = "Cell radius from system centre"
ax3.ylabel = "Aspect ratio"



display(fig)
save(datadir("sims", folderName, "stats.png"), fig)


