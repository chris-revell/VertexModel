using LinearAlgebra
using CairoMakie

#%%

radii = Float64[]
cellAreas = Float64[]
cellAreaErrors = Float64[]
aspectRatios = Float64[]
aspectRatioErrors = Float64[]
hexagonalities = Float64[]
hexagonalityErrors = Float64[]
neighbourCounts = Float64[]
neighbourCountErrors = Float64[]

# P1 area - fovea and parafovea, 3 mm diameter (1.5 mm radius), located at approx. 3 mm from the optic nerve (“blind spot”)
push!(radii, 3.0)
# area 147.24 ± 15.36 µm2, 
push!(cellAreas, 147.24)
push!(cellAreaErrors, 15.36)
# AR 1.15 ± 0.04, 
push!(aspectRatios, 1.15)
push!(aspectRatioErrors, 0.04)
# hexagonality 9.31 ± 0.11, 
push!(hexagonalities, 9.31)
push!(hexagonalityErrors, 0.11)
# number of neighbors 5.55 ± 0.35
push!(neighbourCounts, 5.55)
push!(neighbourCountErrors, 0.35)

# P2 area:  most of the center of the RPE monolayer, including the perifovea, up to a 10-mm radial distance:
# Ave cell area 201.74 ± 17.45 µm2, 
# AR 1.18 ± 0.02, 
# hexagonality 9.25 ± 0.05,
# number of neighbors 5.47 ± 0.43 
push!(radii, 7.5)
push!(cellAreas, 201.74)
push!(cellAreaErrors, 17.45)
push!(aspectRatios, 1.18)
push!(aspectRatioErrors, 0.02)
push!(hexagonalities, 9.25)
push!(hexagonalityErrors, 0.05)
push!(neighbourCounts, 5.47)
push!(neighbourCountErrors, 0.43)

# P3 area: midperipheral ring of RPE cells located 10 to 14 mm from the center
# average cell area of 231.21 ± 18.08 µm2, 
# AR 1.23 ± 0.03, 
# hexagonality 9.12 ± 0.08
# number of neighbors 5.46 ± 0.63
push!(radii, 12.0)
push!(cellAreas, 231.21)
push!(cellAreaErrors, 18.08)
push!(aspectRatios, 1.23)
push!(aspectRatioErrors, 0.03)
push!(hexagonalities, 9.12)
push!(hexagonalityErrors, 0.08)
push!(neighbourCounts, 5.46)
push!(neighbourCountErrors, 0.63)

# P4 area: newly discovered small RPE cells of the periphery situated on average 14 to 17 mm from the center
# average cell area of 176.76 ± 18.68 µm2, 
# AR 1.27 ± 0.04, 
# hexagonality 9.00 ± 0.12, 
# number of neighbors 5.64 ± 0.25
push!(radii, 15.5)
push!(cellAreas, 176.76)
push!(cellAreaErrors, 18.68)
push!(aspectRatios, 1.27)
push!(aspectRatioErrors, 0.04)
push!(hexagonalities, 9.00)
push!(hexagonalityErrors, 0.12)
push!(neighbourCounts, 5.64)
push!(neighbourCountErrors, 0.25)

# P5 area: far-peripheral RPE cells positioned 17 mm away from the center of the flatmount until the start of the ora serrata, at the very edge of the flatmount (22 to 25 mm from the center)
# average cell area of 331.87 ± 27.23 µm2, 
# AR 1.33 ± 0.03, 
# hexagonality 8.79 ± 0.11 
# number of neighbors 5.04 ± 0.46
push!(radii, 17.0)
push!(cellAreas, 331.87)
push!(cellAreaErrors, 27.23)
push!(aspectRatios, 1.33)
push!(aspectRatioErrors, 0.03)
push!(hexagonalities, 8.79)
push!(hexagonalityErrors, 0.11)
push!(neighbourCounts, 5.04)
push!(neighbourCountErrors, 0.46) 


using CairoMakie
allLines = []
axes = []
fig = Figure()
push!(axes, Axis(fig[1,1]))
push!(allLines, lines!(axes[end], radii, cellAreas))
errorbars!(axes[end], radii, cellAreas, cellAreaErrors)
axes[end].ylabel = "cell areas"
axes[end].xlabel = "radius"
push!(axes, Axis(fig[2,1]))
push!(allLines, lines!(axes[end], radii, aspectRatios))
errorbars!(axes[end], radii, aspectRatios, aspectRatioErrors)
axes[end].ylabel = "aspect ratios"
axes[end].xlabel = "radius"
push!(axes, Axis(fig[1,2]))
push!(allLines, lines!(axes[end], radii, hexagonalities))
errorbars!(axes[end], radii, hexagonalities, hexagonalityErrors)
axes[end].ylabel = "hexagonlity"
axes[end].xlabel = "radius"
push!(axes, Axis(fig[2,2]))
push!(allLines, lines!(axes[end], radii, neighbourCounts))
errorbars!(axes[end], radii, neighbourCounts, neighbourCountErrors)
axes[end].ylabel = "neighbour counts"
axes[end].xlabel = "radius"
display(fig)

save("retinaData.png", fig)



