
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

# folderName = "MCCComparison/pressureExternal=0.5_stiffnessFactor=10.0_γ=0.2_24-06-24-21-29-00"

files = [datadir("sims", folderName, "frameData", f) for f in readdir(datadir("sims", folderName, "frameData")) if occursin(".jld2",f)]
@unpack R, matrices, params = load(files[end]; 
    typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer"=>ParametersContainer, 
    "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer"=>MatricesContainer))

series = Dict()

function myECDF!(ax, dataNeighbour, dataBulk, series, label)
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

    series["$(label)MCCneighbours_x"] = xs1
    series["$(label)MCCneighbours_y"] = ys1
    series["$(label)BulkCells_x"] = xs2
    series["$(label)BulkCells_y"] = ys2

    # pts1 = Point2[]
    # push!(pts1, Point2(minVal,0.0))
    # pts2 = Point2[]
    # push!(pts1, Point2(maxVal,0.0))
    # for (i,d) in enumerate(dataNeighbour)
    #     push!(pts1, Point2(d,(i-1)/length(dataNeighbour)))
    #     push!(pts1, Point2(d,i/length(dataNeighbour)))
    # end
    # push!(pts1, Point2((maxVal, 1.0)))
    # for (i,d) in enumerate(dataBulk)
    #     push!(pts2, Point2((d,(i-1)/length(dataBulk))))
    #     push!(pts2, Point2((d,i/length(dataBulk))))
    # end
    # push!(pts2, Point2((maxVal, 1.0)))

    l1 = lines!(ax, xs1, ys1, color=:red, label="MCC neighbours")
    l2 = lines!(ax, xs2, ys2, color=:green, label="Bulk cells")

    return l1, l2
end 

cellNeighbourMatrix = matrices.B*matrices.B'
MCCs = findall(x->x==1, matrices.MCCsList)
MCCneighbours = Int64[]
for i in MCCs
    append!(MCCneighbours, [x for x in findall(x->x!=0, cellNeighbourMatrix[i,:]) if x∉MCCs] )
end
unique!(MCCneighbours)
MCCtensions = matrices.cellTensions[MCCs]
MCCneighbourTensions = matrices.cellTensions[MCCneighbours]
otherTensions = matrices.cellTensions[Not([MCCneighbours...,MCCs...])]
Qs = cellQs(matrices.cellPerimeters, matrices.edgeTangents, matrices.B̄)
shrs = cellShears(matrices.cellTensions, matrices.cellPerimeters, matrices.cellAreas, Qs)
effectivePressure = effectiveCellPressure(matrices.cellPressures, matrices.cellTensions, matrices.cellPerimeters, matrices.cellAreas)

#%%

fig = Figure(size=(1000,1000), fontsize=24); ax = Axis(fig[1,1])
# sc1 = ecdfplot!(ax, matrices.cellTensions[MCCneighbours]; color=:red, label="Neighbour cell tensions")
# sc2 = ecdfplot!(ax, matrices.cellTensions[Not([MCCneighbours...,MCCs...])]; color=:green, label="Other cell tensions")
# axislegend(ax)#, merge = true, unique = true)
# pts1, pts2 = myECDF!(ax, matrices.cellTensions[MCCneighbours], matrices.cellTensions[Not([MCCneighbours...,MCCs...])], "label")
l1, l2 = myECDF!(ax, matrices.cellTensions[MCCneighbours], matrices.cellTensions[Not([MCCneighbours...,MCCs...])], series, "CellTensions")
# lines!(ax, pts1)#, color=:red, label="MCC neighbours")
# lines!(ax, pts2)#, color=:green, label="Bulk cells")
ax.ylabel = "Cumulative distribution"
ax.xlabel = "Tension"
Legend(fig[1, 2], [l1, l2], ["MCC neighbours", "Bulk cells"])
# display(fig)
# axislegend(ax)
save(datadir("sims", folderName, "tensionCumulativeDensities.png"), fig)

#%%

fig = Figure(size=(1000,1000), fontsize=24); ax = Axis(fig[1,1])
# sc1 = ecdfplot!(ax, matrices.cellTensions[MCCneighbours]; color=:red, label="Neighbour cell tensions")
# sc2 = ecdfplot!(ax, matrices.cellTensions[Not([MCCneighbours...,MCCs...])]; color=:green, label="Other cell tensions")
# axislegend(ax)#, merge = true, unique = true)
# pts1, pts2 = myECDF!(ax, matrices.cellTensions[MCCneighbours], matrices.cellTensions[Not([MCCneighbours...,MCCs...])], "label")
l1, l2 = myECDF!(ax, effectivePressure[MCCneighbours], effectivePressure[Not([MCCneighbours...,MCCs...])], series, "CellEffectivePressures")
# lines!(ax, pts1)#, color=:red, label="MCC neighbours")
# lines!(ax, pts2)#, color=:green, label="Bulk cells")
ax.ylabel = "Cumulative distribution"
ax.xlabel = "Effective pressure"
Legend(fig[1, 2], [l1, l2], ["MCC neighbours", "Bulk cells"])
# display(fig)
# axislegend(ax)
save(datadir("sims", folderName, "peffCumulativeDensities.png"), fig)


#%%
#%%

# fig = Figure(size=(1000,1000)); ax = Axis(fig[1,1])
# ecdfplot!(ax, matrices.cellPressures[MCCneighbours], label="Neighbour cell\n pressures", color=:green)
# ecdfplot!(ax, matrices.cellPressures[Not([MCCneighbours...,MCCs...])], label="Other cell\n pressures", color=:blue)
# # axislegend(ax)#, merge = true, unique = true)
# save(datadir("sims", folderName, "pressureCumulativeDensities.png"), fig)



fig = Figure(size=(1000,1000), fontsize=24); ax = Axis(fig[1,1])
# sc1 = ecdfplot!(ax, matrices.cellTensions[MCCneighbours]; color=:red, label="Neighbour cell tensions")
# sc2 = ecdfplot!(ax, matrices.cellTensions[Not([MCCneighbours...,MCCs...])]; color=:green, label="Other cell tensions")
# axislegend(ax)#, merge = true, unique = true)
# pts1, pts2 = myECDF!(ax, matrices.cellTensions[MCCneighbours], matrices.cellTensions[Not([MCCneighbours...,MCCs...])], "label")
l1, l2 = myECDF!(ax, matrices.cellPressures[MCCneighbours], matrices.cellPressures[Not([MCCneighbours...,MCCs...])], series, "CellPressures")
# lines!(ax, pts1)#, color=:red, label="MCC neighbours")
# lines!(ax, pts2)#, color=:green, label="Bulk cells")
ax.ylabel = "Cumulative distribution"
ax.xlabel = "Cell pressure"
Legend(fig[1, 2], [l1, l2], ["MCC neighbours", "Bulk cells"])
# display(fig)
# axislegend(ax)
save(datadir("sims", folderName, "pressureCumulativeDensities.png"), fig)


#%%

# fig = Figure(size=(1000,1000)); ax = Axis(fig[1,1])
# ecdfplot!(ax, shrs[MCCneighbours], label="Neighbour cell\nshears", color=:green)
# ecdfplot!(ax, shrs[Not([MCCneighbours...,MCCs...])], label="Other cell\nshears", color=:blue)
# # axislegend(ax)#, merge = true, unique = true)
# save(datadir("sims", folderName, "shearCumulativeDensities.png"), fig)



fig = Figure(size=(1000,1000), fontsize=24); ax = Axis(fig[1,1])
# sc1 = ecdfplot!(ax, matrices.cellTensions[MCCneighbours]; color=:red, label="Neighbour cell tensions")
# sc2 = ecdfplot!(ax, matrices.cellTensions[Not([MCCneighbours...,MCCs...])]; color=:green, label="Other cell tensions")
# axislegend(ax)#, merge = true, unique = true)
# pts1, pts2 = myECDF!(ax, matrices.cellTensions[MCCneighbours], matrices.cellTensions[Not([MCCneighbours...,MCCs...])], "label")
l1, l2 = myECDF!(ax, shrs[MCCneighbours], shrs[Not([MCCneighbours...,MCCs...])], series, "CellShears")
# lines!(ax, pts1)#, color=:red, label="MCC neighbours")
# lines!(ax, pts2)#, color=:green, label="Bulk cells")
ax.ylabel = "Cumulative distribution"
ax.xlabel = "Cell shear"
Legend(fig[1, 2], [l1, l2], ["MCC neighbours", "Bulk cells"])
# display(fig)
# axislegend(ax)
save(datadir("sims", folderName, "shearCumulativeDensities.png"), fig)


#%%


# using XLSX
# XLSX.openxlsx(datadir("sims", folderName, "ecdfs.xlsx"), mode="w") do xf
#     sheet = xf[1]
#     # XLSX.rename!(sheet, "new_sheet")
#     for (i,j) in enumerate(keys(series))
#         sheet["A$(i)"] = j 
#         sheet["B$(i)"] = series[j]
#     end
# end


# #%%

# cellTypes = fill("Bulk cell", params.nCells)
# cellTypes[MCCneighbours] .= "MCC neighbour"
# dfCells = DataFrame(cellID = collect(1:params.nCells), cellType=cellTypes, pressure=matrices.cellPressures, pEff=effectivePressure, tension=matrices.cellTensions, shear=shrs)
# XLSX.writetable(datadir("sims", folderName, "MCCcellData.xlsx"), collect(eachcol(dfCells)), names(dfCells))