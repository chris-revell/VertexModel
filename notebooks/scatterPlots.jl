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
figTension = Figure(size=(1000,1000), fontsize=24)
axTension = Axis(figTension[1,1], aspect=1)    
allTensions = Float64[]

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
    # allPressures = vcat(N₀PressurePre,N₁PressurePre,N₂PressurePre,N₀PressurePost,N₁PressurePost,N₂PressurePost)
    append!(allPressures, N₀PressurePre,N₁PressurePre,N₂PressurePre,N₀PressurePost,N₁PressurePost,N₂PressurePost)
    scatter!(axPressure, Point{2,Float64}(N₀PressurePre,N₀PressurePost), color=:red, label="N₀", marker=markerstyles[ii])
    scatter!(axPressure, Point{2,Float64}.(N₁PressurePre,N₁PressurePost), color=:green, label="N₁", marker=markerstyles[ii])
    scatter!(axPressure, Point{2,Float64}.(N₂PressurePre,N₂PressurePost), color=:blue, label="N₂", marker=markerstyles[ii])
    axislegend(axPressure)
    axPressure.xlabel = "Cell pressure before size reduction"
    axPressure.ylabel = "Cell pressure after size reduction"
    axPressure.title = "Change in cell pressure after N₀ \ncell preferred area reduction"

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
    axislegend(axPeff)
    axPeff.xlabel = "Cell effective pressure before size reduction"
    axPeff.ylabel = "Cell effective pressure after size reduction"
    axPeff.title = "Change in cell effective pressure after N₀ \ncell preferred area reduction"
    

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

    append!(allTensions,radialEdgeTensionsPre,radialEdgeTensionsPost)
    scatter!(axTension, Point{2,Float64}.(radialEdgeTensionsPre,radialEdgeTensionsPost), color=:red, marker=markerstyles[ii])
    # axislegend(axTension)
    axTension.xlabel = "Radial edge tension before size reduction"
    axTension.ylabel = "Radial edge tension after size reduction"
    axTension.title = "Tensions in radial edges around N₀ before and after area reduction"

end

lines!(axPressure, [Point2(minimum(allPressures),minimum(allPressures)), Point2(maximum(allPressures),maximum(allPressures))], color=(:black,0.5))
lines!(axPeff, [Point2(minimum(allPeffs),minimum(allPeffs)), Point2(maximum(allPeffs),maximum(allPeffs))], color=(:black,0.5))
lines!(axTension, [Point2(minimum(allTensions),minimum(allTensions)), Point2(maximum(allTensions),maximum(allTensions))], color=(:black,0.5))
save(datadir("sims", "cellPressureScatter.png"), figPressure)
save(datadir("sims", "cellPeffScatter.png"), figPeff)
save(datadir("sims", "radialEdgeScatter.png"), figTension)