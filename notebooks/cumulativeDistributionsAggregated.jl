
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

function myECDF!(ax, dataNeighbour, dataBulk, allSeries, label)
    minVal = min(minimum(dataNeighbour), minimum(dataBulk))
    maxVal = max(maximum(dataNeighbour), maximum(dataBulk))

    xs1 = Float64[]
    ys1 = Float64[]
    push!(xs1, minVal)
    push!(ys1, 0.0)
    for (i,d) in enumerate(sort(dataNeighbour))
        push!(xs1, d)
        push!(ys1, (i-1)/length(dataNeighbour))
        push!(xs1, d)
        push!(ys1, i/length(dataNeighbour))
    end
    push!(xs1, maxVal)
    push!(ys1, 1.0)
    
    xs2 = Float64[]
    ys2 = Float64[]
    push!(xs2, minVal)
    push!(ys2, 0.0)
    for (i,d) in enumerate(sort(dataBulk))
        push!(xs2, d)
        push!(ys2, (i-1)/length(dataBulk))
        push!(xs2, d)
        push!(ys2, i/length(dataBulk))
    end
    push!(xs2, maxVal)
    push!(ys2, 1.0)

    allSeries["$(label)MCCneighbours_x"] = xs1
    allSeries["$(label)MCCneighbours_y"] = ys1
    allSeries["$(label)BulkCells_x"] = xs2
    allSeries["$(label)BulkCells_y"] = ys2

    l1 = lines!(ax, xs1, ys1, color=:red, label="MCC neighbours")
    l2 = lines!(ax, xs2, ys2, color=:green, label="Bulk cells")

    return l1, l2
end 

aggregatedBulkTensions = Float64[]
aggregatedMCCneighbourTensions = Float64[]
aggregatedBulkPressures = Float64[]
aggregatedMCCneighbourPressures = Float64[]
aggregatedBulkpEffs = Float64[]
aggregatedMCCneighbourpEffs = Float64[]
aggregatedBulkShears = Float64[]
aggregatedMCCneighbourShears = Float64[]

allSeries = Dict()

for folderName in [r for r in readdir(datadir("sims", "MCCComparison")) if (occursin("stiffnessFactor=10.0_AllDividing",r) && isdir(datadir("sims", "MCCComparison",r)))]

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

end

#%%

fig = Figure(size=(1000,1000), fontsize=24); ax = Axis(fig[1,1])
l1, l2 = myECDF!(ax, aggregatedMCCneighbourTensions, aggregatedBulkTensions, allSeries, "CellTensions")
ax.ylabel = "Cumulative distribution"
ax.xlabel = "Tension"
Legend(fig[1, 2], [l1, l2], ["MCC neighbours", "Bulk cells"])
# axislegend(ax)
save(datadir("sims", "MCCComparison", "tensionCumulativeDensities_stiffnessFactor=10.0_AllDividing.png"), fig)

#%%

fig = Figure(size=(1000,1000), fontsize=24); ax = Axis(fig[1,1])
l1, l2 = myECDF!(ax, aggregatedMCCneighbourpEffs, aggregatedBulkpEffs, allSeries, "CellEffectivePressures")
ax.ylabel = "Cumulative distribution"
ax.xlabel = "Effective pressure"
Legend(fig[1, 2], [l1, l2], ["MCC neighbours", "Bulk cells"])
# axislegend(ax)
save(datadir("sims", "MCCComparison", "peffCumulativeDensities_stiffnessFactor=10.0_AllDividing.png"), fig)


#%%

fig = Figure(size=(1000,1000), fontsize=24); ax = Axis(fig[1,1])
l1, l2 = myECDF!(ax, aggregatedMCCneighbourPressures, aggregatedBulkPressures, allSeries, "CellPressures")
ax.ylabel = "Cumulative distribution"
ax.xlabel = "Cell pressure"
Legend(fig[1, 2], [l1, l2], ["MCC neighbours", "Bulk cells"])
# axislegend(ax)
save(datadir("sims", "MCCComparison", "pressureCumulativeDensities_stiffnessFactor=10.0_AllDividing.png"), fig)

#%%

fig = Figure(size=(1000,1000), fontsize=24); ax = Axis(fig[1,1])
l1, l2 = myECDF!(ax, aggregatedMCCneighbourShears, aggregatedBulkShears, allSeries, "CellShears")
ax.ylabel = "Cumulative distribution"
ax.xlabel = "Cell shear"
Legend(fig[1, 2], [l1, l2], ["MCC neighbours", "Bulk cells"])
# axislegend(ax)
save(datadir("sims", "MCCComparison", "shearCumulativeDensities_stiffnessFactor=10.0_AllDividing.png"), fig)

#%%

letters = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]

using XLSX
XLSX.openxlsx(datadir("sims", "MCCComparison", "ecdfs_stiffnessFactor=10.0_AllDividing.xlsx"), mode="w") do xf
    sheet = xf[1]
    # XLSX.rename!(sheet, "new_sheet")
    for (i,j) in enumerate(keys(allSeries))
        sheet["$(letters[i])1"] = j 
        sheet["$(letters[i])2", dim=1] = allSeries[j]
    end
end


#%%

# simulationLabels = [fill(1, 800); fill(2, 800); fill(3,800)]
# cellTypes = fill("MCC neighbour", length(a))

# cellTypes = fill("Bulk cell", params.nCells)
# cellTypes[MCCneighbours] .= "MCC neighbour"
# dfCells = DataFrame(cellID = collect(1:params.nCells), cellType=cellTypes, pressure=matrices.cellPressures, pEff=effectivePressure, tension=matrices.cellTensions, shear=shrs)
# XLSX.writetable(datadir("sims", "MCCComparison", "MCCcellData_stiffnessFactor=10.0_AllDividing.xlsx"), collect(eachcol(dfCells)), names(dfCells))
