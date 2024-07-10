
using JLD2
using SparseArrays
using StaticArrays
using LinearAlgebra
using DrWatson
using FromFile
using UnPack
using CairoMakie
using Printf
using Colors
using Statistics
using InvertedIndices
using GeometryBasics
using DataFrames

@from "$(projectdir())/src/VertexModelContainers.jl" using VertexModelContainers
@from "$(projectdir())/src/AnalysisFunctions.jl" using AnalysisFunctions

# function myECDF!(data)
#     minVal = minimum(data)
#     maxVal = maximum(data)

#     xs = Float64[]
#     ys = Int64[]
#     push!(xs, minVal)
#     push!(ys, 0)
#     for (i,d) in enumerate(sort(dataNeighbour))
#         push!(xs1, d)
#         push!(ys1, (i-1)/length(dataNeighbour))
#         push!(xs1, d)
#         push!(ys1, i/length(dataNeighbour))
#     end
#     push!(xs1, maxVal)
#     push!(ys1, 1.0)
    
#     xs2 = Float64[]
#     ys2 = Float64[]
#     push!(xs2, minVal)
#     push!(ys2, 0.0)
#     for (i,d) in enumerate(sort(dataBulk))
#         push!(xs2, d)
#         push!(ys2, (i-1)/length(dataBulk))
#         push!(xs2, d)
#         push!(ys2, i/length(dataBulk))
#     end
#     push!(xs2, maxVal)
#     push!(ys2, 1.0)

#     allSeries["$(label)MCCneighbours_x"] = xs1
#     allSeries["$(label)MCCneighbours_y"] = ys1
#     allSeries["$(label)BulkCells_x"] = xs2
#     allSeries["$(label)BulkCells_y"] = ys2

#     l1 = lines!(ax, xs1, ys1, color=:red, label="MCC neighbours")
#     l2 = lines!(ax, xs2, ys2, color=:green, label="Bulk cells")

#     return l1, l2
# end 

allShears = Float64[]
allTensions = Float64[]
allPressures = Float64[]
allEffectivePressures = Float64[]
allSimIDs = Int[]
allCellIDs = Int[]
allCellTypes = String[]
allMCCStiffnessFactors = String[]
allMCCDivisionStates = String[]

setupparams = [("1.0", "AllDividing"), 
    ("10.0", "AllDividing"),
    ("1.0", "MCCs not dividing"),
    ("10.0", "MCCs not dividing"),
    ("2.0", "MCCs not dividing"),
]

for p in setupparams

    stiffnessFactor, divisionState = p
    
    if divisionState=="AllDividing"
        occurString = "stiffnessFactor=$(stiffnessFactor)_$(divisionState)_γ"
    else
        occurString = "stiffnessFactor=$(stiffnessFactor)_γ"
    end
        
    allSeries = Dict()
    aggregatedBulkTensions = Float64[]
    aggregatedMCCneighbourTensions = Float64[]
    aggregatedBulkPressures = Float64[]
    aggregatedMCCneighbourPressures = Float64[]
    aggregatedBulkpEffs = Float64[]
    aggregatedMCCneighbourpEffs = Float64[]
    aggregatedBulkShears = Float64[]
    aggregatedMCCneighbourShears = Float64[]

    figAllTensionsLines = Figure(size=(1000,1000), fontsize=24)
    axAllTensionsLines = Axis(figAllTensionsLines[1,1])
    figAllpEffsLines = Figure(size=(1000,1000), fontsize=24)
    axAllpEffsLines = Axis(figAllpEffsLines[1,1])
    figAllPressuresLines = Figure(size=(1000,1000), fontsize=24)
    axAllPressuresLines = Axis(figAllPressuresLines[1,1])
    figAllShearLines = Figure(size=(1000,1000), fontsize=24)
    axAllShearLines = Axis(figAllShearLines[1,1])

    for folderName in [r for r in readdir(datadir("sims", "MCCComparison")) if (occursin(occurString,r) && isdir(datadir("sims", "MCCComparison",r)))]

        files = [datadir("sims", "MCCComparison", folderName, "frameData", f) for f in readdir(datadir("sims", "MCCComparison", folderName, "frameData")) if occursin(".jld2",f)]
        @unpack R, matrices, params = load(files[end]; 
            typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer"=>ParametersContainer, 
            "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer"=>MatricesContainer))

        peripheralCells = [findall(z->z!=0, matrices.B[:,x])[1] for x in findall(y->y!=0, matrices.boundaryEdges)]

        cellNeighbourMatrix = matrices.B*matrices.B'
        MCCs = findall(x->x==1, matrices.MCCsList)
        MCCneighbours = Int64[]
        for i in MCCs
            append!(MCCneighbours, [x for x in findall(x->x!=0, cellNeighbourMatrix[i,:]) if x∉MCCs] )
        end
        setdiff!(MCCneighbours, peripheralCells)
        unique!(MCCneighbours)
        
        Qs = cellQs(matrices.cellPerimeters, matrices.edgeTangents, matrices.B̄)
        shrs = cellShears(matrices.cellTensions, matrices.cellPerimeters, matrices.cellAreas, Qs)
        effectivePressure = effectiveCellPressure(matrices.cellPressures, matrices.cellTensions, matrices.cellPerimeters, matrices.cellAreas)

        append!(aggregatedBulkTensions, matrices.cellTensions[Not([MCCneighbours...,MCCs...,peripheralCells...])])
        append!(aggregatedMCCneighbourTensions, matrices.cellTensions[MCCneighbours])
        append!(aggregatedBulkPressures, matrices.cellPressures[Not([MCCneighbours...,MCCs...,peripheralCells...])])
        append!(aggregatedMCCneighbourPressures, matrices.cellPressures[MCCneighbours])
        append!(aggregatedBulkpEffs, effectivePressure[Not([MCCneighbours...,MCCs...,peripheralCells...])])
        append!(aggregatedMCCneighbourpEffs, effectivePressure[MCCneighbours])
        append!(aggregatedBulkShears, shrs[Not([MCCneighbours...,MCCs...,peripheralCells...])])
        append!(aggregatedMCCneighbourShears, shrs[MCCneighbours])

        append!(allShears, shrs)
        append!(allTensions, matrices.cellTensions)
        append!(allPressures, matrices.cellPressures)
        append!(allEffectivePressures, effectivePressure)
        cellTypes = fill("Bulk", params.nCells)
        for i in MCCs
            cellTypes[i] = "MCC"
        end
        for i in MCCneighbours
            cellTypes[i] = "MCC neighbour"
        end
        for i in peripheralCells
            cellTypes[i] = "Peripheral cell"
        end
        append!(allCellTypes, cellTypes)
        append!(allMCCStiffnessFactors, fill(stiffnessFactor, params.nCells))
        append!(allMCCDivisionStates, fill(divisionState, params.nCells))

        xs = sort(matrices.cellTensions[Not([MCCneighbours...,MCCs...,peripheralCells...])])
        ys = collect(1:length(xs))./length(xs)
        stairs!(axAllTensionsLines, xs, ys, label="Bulk Tensions", color=(:green,0.5))
        xs = sort(matrices.cellTensions[MCCneighbours])
        ys = collect(1:length(xs))./length(xs)
        stairs!(axAllTensionsLines, xs, ys, label="Neighbour Tensions", color=(:red,0.5))
    
        xs = sort(effectivePressure[Not([MCCneighbours...,MCCs...,peripheralCells...])])
        ys = collect(1:length(xs))./length(xs)
        stairs!(axAllpEffsLines, xs, ys, label="Bulk pEffs", color=(:green,0.5))
        xs = sort(effectivePressure[MCCneighbours])
        ys = collect(1:length(xs))./length(xs)
        stairs!(axAllpEffsLines, xs, ys, label="Neighbour pEffs", color=(:red,0.5))
    
        xs = sort(matrices.cellPressures[Not([MCCneighbours...,MCCs...,peripheralCells...])])
        ys = collect(1:length(xs))./length(xs)
        stairs!(axAllPressuresLines, xs, ys, label="Bulk pressures", color=(:green,0.5))
        xs = sort(matrices.cellPressures[MCCneighbours])
        ys = collect(1:length(xs))./length(xs)
        stairs!(axAllPressuresLines, xs, ys, label="Neighbour pressures", color=(:red,0.5))
        
        xs = sort(shrs[Not([MCCneighbours...,MCCs...,peripheralCells...])])
        ys = collect(1:length(xs))./length(xs)
        stairs!(axAllShearLines, xs, ys, label="Bulk pEffs", color=(:green,0.5))
        xs = sort(shrs[MCCneighbours])
        ys = collect(Float64, 1:length(xs))./length(xs)
        stairs!(axAllShearLines, xs, ys, label="Neighbour pEffs", color=(:red,0.5))

    end

    save(datadir("sims", "MCCComparison", "tensionCumulativeDensities_$(occurString[1:end-2])_allRuns.png"), figAllTensionsLines)
    save(datadir("sims", "MCCComparison", "peffCumulativeDensities_$(occurString[1:end-2])_allRuns.png"), figAllpEffsLines)
    save(datadir("sims", "MCCComparison", "pressureCumulativeDensities_$(occurString[1:end-2])_allRuns.png"), figAllPressuresLines)
    save(datadir("sims", "MCCComparison", "shearCumulativeDensities_$(occurString[1:end-2])_allRuns.png"), figAllShearLines)
    
    #%%

    fig = Figure(size=(1000,1000), fontsize=24); ax = Axis(fig[1,1])
    xs = sort(aggregatedMCCneighbourTensions)
    ys = collect(1:length(xs))./length(xs)
    l1 = stairs!(ax, xs, ys, label="Neighbour Tensions", color=(:red,0.5))
    xs = sort(aggregatedBulkTensions)
    ys = collect(1:length(xs))./length(xs)
    l2 = stairs!(ax, xs, ys, label="Bulk Tensions", color=(:green,0.5))
    ax.ylabel = "Cumulative distribution"
    ax.xlabel = "Tension"
    Legend(fig[1, 2], [l1, l2], ["MCC neighbours", "Bulk cells"])
    save(datadir("sims", "MCCComparison", "tensionCumulativeDensities_$(occurString[1:end-2]).png"), fig)

    #%%

    fig = Figure(size=(1000,1000), fontsize=24); ax = Axis(fig[1,1])
    xs = sort(aggregatedMCCneighbourpEffs)
    ys = collect(1:length(xs))./length(xs)
    l1 = stairs!(ax, xs, ys, label="Neighbour pEffs", color=(:red,0.5))
    xs = sort(aggregatedBulkpEffs)
    ys = collect(1:length(xs))./length(xs)
    l2 = stairs!(ax, xs, ys, label="Bulk pEffs", color=(:green,0.5))
    ax.ylabel = "Cumulative distribution"
    ax.xlabel = "Effective pressure"
    Legend(fig[1, 2], [l1, l2], ["MCC neighbours", "Bulk cells"])
    save(datadir("sims", "MCCComparison", "peffCumulativeDensities_$(occurString[1:end-2]).png"), fig)


    #%%

    fig = Figure(size=(1000,1000), fontsize=24); ax = Axis(fig[1,1])
    xs = sort(aggregatedMCCneighbourPressures)
    ys = collect(1:length(xs))./length(xs)
    l1 = stairs!(ax, xs, ys, label="Neighbour Pressures", color=(:red,0.5))
    xs = sort(aggregatedBulkPressures)
    ys = collect(1:length(xs))./length(xs)
    l2 = stairs!(ax, xs, ys, label="Bulk Pressures", color=(:green,0.5))
    ax.ylabel = "Cumulative distribution"
    ax.xlabel = "Cell pressure"
    Legend(fig[1, 2], [l1, l2], ["MCC neighbours", "Bulk cells"])
    # axislegend(ax)
    save(datadir("sims", "MCCComparison", "pressureCumulativeDensities_$(occurString[1:end-2]).png"), fig)

    #%%

    fig = Figure(size=(1000,1000), fontsize=24); ax = Axis(fig[1,1])
    xs = sort(aggregatedMCCneighbourShears)
    ys = collect(1:length(xs))./length(xs)
    l1 = stairs!(ax, xs, ys, label="Neighbour Shears", color=(:red,0.5))
    xs = sort(aggregatedBulkShears)
    ys = collect(1:length(xs))./length(xs)
    l2 = stairs!(ax, xs, ys, label="Bulk Shears", color=(:green,0.5))
    ax.ylabel = "Cumulative distribution"
    ax.xlabel = "Cell shear"
    Legend(fig[1, 2], [l1, l2], ["MCC neighbours", "Bulk cells"])
    # axislegend(ax)
    save(datadir("sims", "MCCComparison", "shearCumulativeDensities_$(occurString[1:end-2]).png"), fig)

    #%%

    # letters = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
    # using XLSX
    # XLSX.openxlsx(datadir("sims", "MCCComparison", "ecdfs_$(occurString[1:end-2]).xlsx"), mode="w") do xf
    #     sheet = xf[1]
    #     # XLSX.rename!(sheet, "new_sheet")
    #     for (i,j) in enumerate(keys(allSeries))
    #         sheet["$(letters[i])1"] = j 
    #         sheet["$(letters[i])2", dim=1] = allSeries[j]
    #     end
    # end


end
#%%

# simulationLabels = [fill(1, 800); fill(2, 800); fill(3,800)]
# cellTypes = fill("MCC neighbour", length(a))

# cellTypes = fill("Bulk cell", params.nCells)
# cellTypes[MCCneighbours] .= "MCC neighbour"
# df = DataFrame(Shears = allShears,
#     Tensions = allTensions,
#     Pressures = allPressures,
#     EffectivePressures = allEffectivePressures,
#     # CellIDs = allCellIDs,
#     CellTypes = allCellTypes,
#     MCCStiffnessFactors = allMCCStiffnessFactors,
#     MCCDivisionStates = allMCCDivisionStates,
# )

# XLSX.writetable(datadir("sims", "MCCComparison", "allRawData.xlsx"), collect(eachcol(df)), names(df))
