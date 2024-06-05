# Script to produce a movie of Airy Stress evolution over cell network for a given system

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
using DataFrames
using XLSX

@from "$(projectdir())/src/VertexModelContainers.jl" using VertexModelContainers
@from "$(projectdir())/src/AnalysisFunctions.jl" using AnalysisFunctions

# folderNames = ["L₀=0.75_nCells=61_pressureExternal=0.5_realTimetMax=86400.0_stiffnessFactor=1.0_γ=0.2_24-05-24-14-37-24",
#                "L₀=0.75_nCells=61_pressureExternal=0.5_realTimetMax=86400.0_stiffnessFactor=1.0_γ=0.2_24-05-08-16-38-39",
#                "L₀=0.75_nCells=61_pressureExternal=0.5_realTimetMax=86400.0_stiffnessFactor=1.0_γ=0.2_24-05-24-14-35-43",
#                "L₀=0.75_nCells=61_pressureExternal=0.5_realTimetMax=86400.0_stiffnessFactor=1.0_γ=0.2_24-05-24-17-24-44",
#                "L₀=0.75_nCells=61_pressureExternal=0.5_realTimetMax=86400.0_stiffnessFactor=1.0_γ=0.2_24-05-24-17-25-25"]

folderNames = [f for f in readdir(datadir("sims")) if occursin("L₀",f)]

markerstyles = [:circle,
    :rect,
    :diamond,
    :hexagon,
    :cross,
    :xcross,
    :utriangle,
    :dtriangle,
    :ltriangle,
    :rtriangle,
    :pentagon]

N₀ = 5

figPressure = Figure(size=(1000,1000), fontsize=24)
axPressure = Axis(figPressure[1,1], aspect=1)
allPressures = Float64[]
figPeff = Figure(size=(1000,1000), fontsize=24)
axPeff = Axis(figPeff[1,1], aspect=1)
allPeffs = Float64[]
figCellTension = Figure(size=(1000,1000), fontsize=24)
axCellTension = Axis(figCellTension[1,1], aspect=1)    
allCellTensions = Float64[]
figEdgeTension = Figure(size=(1000,1000), fontsize=24)
axEdgeTension = Axis(figEdgeTension[1,1], aspect=1)    
allEdgeTensions = Float64[]


dfCells = DataFrame(cellID = Int64[], simulationID = Int64[], cellType=String[], pressurePre=Float64[], pressurePost=Float64[], pEffPre=Float64[], pEffPost=Float64[], tensionPre=Float64[], tensionPost=Float64[])
dfEdges = DataFrame(edgeID = Int64[], simulationID = Int64[], tensionPre=Float64[], tensionPost=Float64[])

# cellID = [1, 2, 3, 4]
# simulationID = [5, 6, 7, 8]
# cellType=["N_0", "N_0","N_0", "N_1"]
# pressurePre=[0.1, 0.2, 0.3, 0.4]
# pressurePost=[0.5, 0.6, 0.7, 0.8]
# pEffPre=[0.1, 0.2, 0.3, 0.4]
# pEffPost=[0.5, 0.6, 0.7, 0.8]
# tensionPre=[0.1, 0.2, 0.3, 0.4]
# tensionPost=[0.5, 0.6, 0.7, 0.8]

# test = Dict(names(dfCells).=>(cellID,simulationID,cellType,pressurePre,pressurePost,pEffPre,pEffPost,tensionPre,tensionPost))

# append!(dfCells, test)
# append!(dfCells, Dict(names(dfCells).=>(cellID,simulationID,cellType,pressurePre,pressurePost,pEffPre,pEffPost,tensionPre,tensionPost)))


for (ii,folderName) in enumerate(folderNames)

    files = [datadir("sims", folderName, "frameData", f) for f in readdir(datadir("sims", folderName, "frameData")) if occursin(".jld2",f)]

    dataPre = load(files[70]; 
            typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer"=>ParametersContainer, 
            "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer"=>MatricesContainer))

    dataPost = load(files[end]; 
            typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer"=>ParametersContainer, 
            "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer"=>MatricesContainer))

    cellNeighbourMatrixPre = dataPre["matrices"].B*dataPre["matrices"].B'
    N₁Pre = [x for x in findall(x->x!=0, cellNeighbourMatrixPre[N₀,:]) if x!=N₀]
    cellNeighbourMatrixPost = dataPost["matrices"].B*dataPost["matrices"].B'
    N₁Post = [x for x in findall(x->x!=0, cellNeighbourMatrixPost[N₀,:]) if x!=N₀]

    N₁Pre!=N₁Post ? println("$folderName") : nothing

    N₂Pre = Int64[]
    for i in N₁Pre 
        iNeighbours = [x for x in findall(x->x!=0, cellNeighbourMatrixPre[i,:]) if (x∉N₁Pre && x!=N₀)]
        append!(N₂Pre, iNeighbours)
    end
    unique!(N₂Pre)
    N₂Post = Int64[]
    for i in N₁Post 
        iNeighbours = [x for x in findall(x->x!=0, cellNeighbourMatrixPost[i,:]) if (x∉N₁Post && x!=N₀)]
        append!(N₂Post, iNeighbours)
    end
    unique!(N₂Post)

    N₂Pre!=N₂Post ? println("$folderName") : nothing

    N₁Both = N₁Pre∩N₁Post
    N₂Both = N₂Pre∩N₂Post
    
    N₀PressurePre = dataPre["matrices"].cellPressures[N₀]
    N₁PressurePre = dataPre["matrices"].cellPressures[N₁Both]
    N₂PressurePre = dataPre["matrices"].cellPressures[N₂Both]
    N₀PressurePost = dataPost["matrices"].cellPressures[N₀]
    N₁PressurePost = dataPost["matrices"].cellPressures[N₁Both]
    N₂PressurePost = dataPost["matrices"].cellPressures[N₂Both]
    append!(allPressures, N₀PressurePre,N₁PressurePre,N₂PressurePre,N₀PressurePost,N₁PressurePost,N₂PressurePost)
    scatter!(axPressure, Point{2,Float64}(N₀PressurePre,N₀PressurePost), color=:red, label="N₀", marker=markerstyles[ii])
    scatter!(axPressure, Point{2,Float64}.(N₁PressurePre,N₁PressurePost), color=:green, label="N₁", marker=markerstyles[ii])
    scatter!(axPressure, Point{2,Float64}.(N₂PressurePre,N₂PressurePost), color=:blue, label="N₂", marker=markerstyles[ii])
    
    N₀TensionPre = dataPre["matrices"].cellTensions[N₀]
    N₁TensionPre = dataPre["matrices"].cellTensions[N₁Both]
    N₂TensionPre = dataPre["matrices"].cellTensions[N₂Both]
    N₀TensionPost = dataPost["matrices"].cellTensions[N₀]
    N₁TensionPost = dataPost["matrices"].cellTensions[N₁Both]
    N₂TensionPost = dataPost["matrices"].cellTensions[N₂Both]
    append!(allCellTensions, N₀TensionPre,N₁TensionPre,N₂TensionPre,N₀TensionPost,N₁TensionPost,N₂TensionPost)
    scatter!(axCellTension, Point{2,Float64}(N₀TensionPre,N₀TensionPost), color=:red, label="N₀", marker=markerstyles[ii])
    scatter!(axCellTension, Point{2,Float64}.(N₁TensionPre,N₁TensionPost), color=:green, label="N₁", marker=markerstyles[ii])
    scatter!(axCellTension, Point{2,Float64}.(N₂TensionPre,N₂TensionPost), color=:blue, label="N₂", marker=markerstyles[ii])
    
    effectivePressurePre = effectiveCellPressure(dataPre["matrices"].cellPressures, dataPre["matrices"].cellTensions, dataPre["matrices"].cellPerimeters, dataPre["matrices"].cellAreas)
    effectivePressurePost = effectiveCellPressure(dataPost["matrices"].cellPressures, dataPost["matrices"].cellTensions, dataPost["matrices"].cellPerimeters, dataPost["matrices"].cellAreas)
    N₀PeffPre = effectivePressurePre[N₀]
    N₁PeffPre = effectivePressurePre[N₁Both]
    N₂PeffPre = effectivePressurePre[N₂Both]
    N₀PeffPost = effectivePressurePost[N₀]
    N₁PeffPost = effectivePressurePost[N₁Both]
    N₂PeffPost = effectivePressurePost[N₂Both]
    append!(allPeffs,N₀PeffPre,N₁PeffPre,N₂PeffPre,N₀PeffPost,N₁PeffPost,N₂PeffPost)
    scatter!(axPeff, Point{2,Float64}(N₀PeffPre,N₀PeffPost), color=:red, label="N₀", marker=markerstyles[ii])
    scatter!(axPeff, Point{2,Float64}.(N₁PeffPre,N₁PeffPost), color=:green, label="N₁", marker=markerstyles[ii])
    scatter!(axPeff, Point{2,Float64}.(N₂PeffPre,N₂PeffPost), color=:blue, label="N₂", marker=markerstyles[ii])
    
    append!(dfCells, Dict(names(dfCells).=>(N₀,ii,"N_0",N₀PressurePre,N₀PressurePost,N₀PeffPre,N₀PeffPost,N₀TensionPre,N₀TensionPost)))
    append!(dfCells, Dict(names(dfCells).=>(N₁Both,fill(ii,length(N₁Both)),fill("N_1",length(N₁Both)),N₁PressurePre,N₁PressurePost,N₁PeffPre,N₁PeffPost,N₁TensionPre,N₁TensionPost)))
    append!(dfCells, Dict(names(dfCells).=>(N₂Both,fill(ii,length(N₂Both)),fill("N_2",length(N₂Both)),N₂PressurePre,N₂PressurePost,N₂PeffPre,N₂PeffPost,N₂TensionPre,N₂TensionPost)))

    N₀VerticesPre = findall(x->x!=0, dataPre["matrices"].C[N₀,:])
    N₀VerticesPost = findall(x->x!=0, dataPre["matrices"].C[N₀,:])
    N₀VerticesPreAndPost = N₀VerticesPre ∩ N₀VerticesPost
    radialEdgesPreAndPost = Int64[]
    for k in N₀VerticesPreAndPost
        radialEdge = [x for x in findall(x->x!=0, dataPre["matrices"].A[:,k]) if dataPre["matrices"].B[N₀,x]==0][1]
        push!(radialEdgesPreAndPost, radialEdge)
    end
    radialEdgeTensionsPre = (dataPre["matrices"].B̄ᵀ*dataPre["matrices"].cellTensions)[radialEdgesPreAndPost]
    radialEdgeTensionsPost = (dataPost["matrices"].B̄ᵀ*dataPost["matrices"].cellTensions)[radialEdgesPreAndPost]

    append!(allEdgeTensions,radialEdgeTensionsPre,radialEdgeTensionsPost)
    scatter!(axEdgeTension, Point{2,Float64}.(radialEdgeTensionsPre,radialEdgeTensionsPost), color=:red, marker=markerstyles[ii])
    
    append!(dfEdges, Dict(names(dfEdges).=>(N₀VerticesPreAndPost,fill(ii,length(N₀VerticesPreAndPost)),radialEdgeTensionsPre,radialEdgeTensionsPost)))
end

#%%

lines!(axPressure, [Point2(minimum(allPressures),minimum(allPressures)), Point2(maximum(allPressures),maximum(allPressures))], color=(:black,0.5))
lines!(axPeff, [Point2(minimum(allPeffs),minimum(allPeffs)), Point2(maximum(allPeffs),maximum(allPeffs))], color=(:black,0.5))
lines!(axEdgeTension, [Point2(minimum(allEdgeTensions),minimum(allEdgeTensions)), Point2(maximum(allEdgeTensions),maximum(allEdgeTensions))], color=(:black,0.5))
lines!(axCellTension, [Point2(minimum(allCellTensions),minimum(allCellTensions)), Point2(maximum(allCellTensions),maximum(allCellTensions))], color=(:black,0.5))
axPressure.xlabel = "Cell pressure before size reduction"
axPressure.ylabel = "Cell pressure after size reduction"
axPressure.title = "Change in cell pressure after N₀ \ncell preferred area reduction"
axCellTension.xlabel = "Cell Tension before size reduction"
axCellTension.ylabel = "Cell Tension after size reduction"
axCellTension.title = "Change in cell tensions after N₀ \ncell preferred area reduction"
axPeff.xlabel = "Cell effective pressure before size reduction"
axPeff.ylabel = "Cell effective pressure after size reduction"
axPeff.title = "Change in cell effective pressure after N₀ \ncell preferred area reduction"
axEdgeTension.xlabel = "Radial edge tension before size reduction"
axEdgeTension.ylabel = "Radial edge tension after size reduction"
axEdgeTension.title = "Tensions in radial edges around N₀ before and after area reduction"
save(datadir("sims", "cellPressureScatter.png"), figPressure)
save(datadir("sims", "cellPeffScatter.png"), figPeff)
save(datadir("sims", "radialEdgeScatter.png"), figEdgeTension)
save(datadir("sims", "cellTenionScatter.png"), figCellTension)


#%%

XLSX.writetable(datadir("sims", "cellData.xlsx"), collect(eachcol(dfCells)), names(dfCells))
XLSX.writetable(datadir("sims", "edgeData.xlsx"), collect(eachcol(dfEdges)), names(dfEdges))