# Script to plot histograms of the recoil speed for each edge type. 

using DrWatson
using DiscreteCalculus
using CairoMakie
using StaticArrays
using VertexModel
using LinearAlgebra
using SparseArrays
using Random
using Colors 
using JLD2
using Dates
using FromFile
using InvertedIndices
using LaTeXStrings
using Dates 
using CircularArrays
using OrdinaryDiffEq
using Printf
using DiffEqCallbacks


# Pick the equilibrated system: 
dateString = "26-09-10-10-49-26"

!isdir(datadir("multipleRuns", dateString)) ? mkpath(datadir("multipleRuns", dateString)) : nothing 

parameterLabelVec = ["(O)","(P)","(Q)","(R)","(S)",
                    "(T)",
                    "(U)","(V)"]

jld2pathVec = ["data/multipleRuns/26-09-10-10-49-26/(O)/26-09-12-14-42-03_nCells=1273_Λ_AA=-0.2_Λ_AB=-0.2_Λ_BB=-0.2_β=0.0_γ=0.05/frameData/systemData099.jld2", #(O)
                "data/multipleRuns/26-09-10-10-49-26/(P)/26-09-12-13-53-32_nCells=1273_Λ_AA=-0.1_Λ_AB=-0.2_Λ_BB=-0.25_β=0.0_γ=0.05/frameData/systemData099.jld2",#(P)
                "data/multipleRuns/26-09-10-10-49-26/(Q)/26-09-12-18-25-13_nCells=1273_Λ_AA=-0.1_Λ_AB=-0.2_Λ_BB=-0.3_β=0.0_γ=0.05/frameData/systemData099.jld2",#(Q)
                "data/multipleRuns/26-09-10-10-49-26/(R)/26-09-12-13-39-48_nCells=1273_Λ_AA=-0.15_Λ_AB=-0.2_Λ_BB=-0.3_β=0.0_γ=0.05/frameData/systemData099.jld2",#(R)
                "data/multipleRuns/26-09-10-10-49-26/(S)/26-09-13-03-46-48_nCells=1169_Λ_AA=-0.3_Λ_AB=-0.2_Λ_BB=-0.3_β=0.0_γ=0.05/frameData/systemData099.jld2",#(S)
                "data/multipleRuns/26-09-10-10-49-26/(T)/26-09-12-14-02-50_nCells=1273_Λ_AA=-0.3_Λ_AB=-0.2_Λ_BB=-0.15_β=0.0_γ=0.05/frameData/systemData099.jld2",#(T)
                "data/multipleRuns/26-09-10-10-49-26/(U)/26-09-12-14-08-47_nCells=1273_Λ_AA=-0.3_Λ_AB=-0.2_Λ_BB=-0.1_β=0.0_γ=0.05/frameData/systemData099.jld2",#(U)
                "data/multipleRuns/26-09-10-10-49-26/(V)/26-09-12-13-57-33_nCells=1273_Λ_AA=-0.25_Λ_AB=-0.2_Λ_BB=-0.1_β=0.0_γ=0.05/frameData/systemData099.jld2",]#(V)


set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
fig = Figure(size=(1800,600))

# Initialise figure: 
grid = fig[1,1] = GridLayout()
ax = Axis(grid[1,1],
xticks = (1:length(parameterLabelVec),parameterLabelVec),
xlabel = "Simulation label",
ylabel = "Recoil speed")

AAcolour = RGB(255/255, 178/255, 102/255)
BBcolour = RGB(102/255, 178/255, 255/255)

replotWithNewColourBar = true


for (simInd,parameterLabel) in enumerate(parameterLabelVec)

    matrixData = jld2pathVec[simInd]
    matrices = load(matrixData,"matrices")
    params = load(matrixData,"params")
    R = load(matrixData,"R")

    jld2pathString = datadir("multipleRuns",dateString,parameterLabel,"ablationLoop","ablationData.jld2")
    edgeSpeeds = load(jld2pathString,"edgeSpeeds")
    edgeRotations = load(jld2pathString,"edgeRotations")

    @unpack edgeLabels = matrices

    typeAAspeeds = []
    typeBBspeeds = []
    typeABspeeds = []

    for j in eachindex(edgeLabels)
        if edgeLabels[j] ==0
            push!(typeAAspeeds,edgeSpeeds[j])
        elseif edgeLabels[j]==1
            push!(typeBBspeeds,edgeSpeeds[j])
        elseif edgeLabels[j]==2
            push!(typeABspeeds,edgeSpeeds[j])
        end
    end

    boxplot!(ax, fill(simInd, length(typeAAspeeds)), typeAAspeeds,
            dodge = fill(1, length(typeAAspeeds)), n_dodge = 3, color =AAcolour, width = 0.7)
    boxplot!(ax, fill(simInd, length(typeBBspeeds)), typeBBspeeds,
            dodge = fill(2, length(typeBBspeeds)), n_dodge = 3, color =BBcolour, width = 0.7)
    boxplot!(ax, fill(simInd, length(typeABspeeds)), typeABspeeds,
            dodge = fill(3, length(typeABspeeds)), n_dodge = 3, color =:grey, width = 0.7)

    if replotWithNewColourBar

        @unpack nEdges = params 
        @unpack A, boundaryEdges = matrices

        set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
        edgeAblationFig = Figure(size=(1200,600))
        grid = edgeAblationFig[1,1] = GridLayout()
        rotationAx = Axis(grid[1,1],aspect=DataAspect())
        speedAx = Axis(grid[1,2],aspect=DataAspect())

        hidedecorations!(rotationAx)
        hidedecorations!(speedAx)
        hidespines!(rotationAx)
        hidespines!(speedAx)

        cmapRotation = cgrad([
                RGB(0.0, 0.0, 1.0),    # blue
                RGB(1.0, 1.0, 1.0),   # white, zero
                RGB(1.0, 0.0, 0.0)   # red
            ], 256)
        climsRotation = (-maximum(abs.(edgeRotations)), maximum(abs.(edgeRotations)))
        # Colour bar exclusing exterior vertices 

        cmapSpeed = cgrad([
                :magenta,
                RGB(0.15, 0.15, 0.15),
                :cyan
            ])
        # climsSpeed = (minimum(edgeSpeeds), maximum(edgeSpeeds))
        nonzeroSpeeds = filter(!iszero, edgeSpeeds)
        climsSpeed = (minimum(nonzeroSpeeds), maximum(nonzeroSpeeds))

        nanPt = Point2f(NaN, NaN)
        edgeVectors = fill((nanPt, nanPt), nEdges)   # length nEdges, all "empty" segments

        for j in 1:nEdges
            boundaryEdges[j] == 1 && continue         # stays NaN, so it isn't drawn
            verts = findall(!=(0), @view A[j, :])
            edgeVectors[j] = (Point2f(R[verts[1]]), Point2f(R[verts[2]]))
        end

        pts = [p for seg in edgeVectors for p in seg if !isnan(p[1])]
        xs = [p[1] for p in pts]
        ys = [p[2] for p in pts]

        xlims!(rotationAx, minimum(xs), maximum(xs))
        ylims!(rotationAx, minimum(ys), maximum(ys))
        xlims!(speedAx, minimum(xs), maximum(xs))
        ylims!(speedAx, minimum(ys), maximum(ys))


        linesegments!(rotationAx, edgeVectors; color = edgeRotations,colorrange = climsRotation, colormap = cmapRotation, linewidth=3)
        linesegments!(speedAx, edgeVectors; color = edgeSpeeds,colorrange = climsSpeed, colormap = cmapSpeed, linewidth=3)

        colourBarFig = Figure(size=(1200,600))
        colourBarGrid = colourBarFig[1,1] = GridLayout()

        cbarRotation = Colorbar(colourBarGrid[2,1],colormap = cmapRotation, colorrange=climsRotation, label="Edge rotation", width=20,height=Relative(0.6))
        cbarSpeed = Colorbar(colourBarGrid[2,2],colormap = cmapSpeed, colorrange=climsSpeed, label="Edge Speed", width=20,height=Relative(0.6))

        println("reaches here on round $parameterLabel")

        save(datadir("multipleRuns", dateString, parameterLabel,"ablationLoop", "colourBar.png"),colourBarFig)
        save(datadir("multipleRuns", dateString, parameterLabel,"ablationLoop", "ablationFigure.png"),edgeAblationFig)
    end

end

save(datadir("multipleRuns", dateString,"recoilSpeedBoxPlot.png"),fig)