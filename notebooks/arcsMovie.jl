
using JLD2
using SparseArrays
using StaticArrays
using LinearAlgebra
using DrWatson
using FromFile
using UnPack
using Printf
using ColorSchemes
using Colors
using GeometryBasics
using Random
using CairoMakie

@from "$(projectdir())/src/VertexModelContainers.jl" using VertexModelContainers
@from "$(projectdir())/src/AnalysisFunctions.jl" using AnalysisFunctions

folderName = "L₀=0.75_nCells=61_pressureExternal=0.5_realTimetMax=86400.0_stiffnessFactor=1.0_γ=0.2_24-05-08-16-38-39"

# set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
fig = Figure(size=(1000,1000))
grid = fig[1,1] = GridLayout()
ax1 = Axis(grid[1,1],aspect=DataAspect())
hidedecorations!(ax1)
hidespines!(ax1)
# Create animation object for visualisation
mov = VideoStream(fig, framerate=5)

files = [datadir("sims", folderName, "frameData", f) for f in readdir(datadir("sims", folderName, "frameData")) if occursin(".jld2",f)]    

for i=1:length(files) 

    @unpack R, matrices, params = load(files[i]; 
    typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer"=>ParametersContainer, 
    "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer"=>MatricesContainer))

    @unpack boundaryEdges,
        A,
        Ā,
        B,
        B̄,
        cellPositions,
        cellPressures,
        cellTensions,
        edgeTangents,
        edgeLengths,
        edgeMidpoints,
        ϵ,
        boundaryVertices = matrices
    @unpack nEdges,
        nVerts,
        nCells = params


    empty!(ax1)

    ax1.title = "t = $(@sprintf("%.2f", (i-1)*params.outputInterval))"

    for j=1:params.nEdges 
        if boundaryEdges[j] == 0
            i,ii = findall(x->x!=0,B[:,j])
            t̲ = edgeTangents[j]
            t = norm(t̲)
            t̂ = t̲./t
            k = findall(x->x<0, A[j,:])
            resultantF = 0.5*cellPressures[i]*B[i,j]*Ā[j,k].*(ϵ*edgeTangents[j])
            resultantF += 0.5*cellPressures[ii]*B[ii,j]*Ā[j,k].*(ϵ*edgeTangents[j])
            resultantF -= cellTensions[i]*B̄[i,j]*A[j,k].*edgeTangents[j]./edgeLengths[j]
            resultantF -= cellTensions[ii]*B̄[ii,j]*A[j,k].*edgeTangents[j]./edgeLengths[j]
            θ = acos(normalize(resultantF)⋅t̂)
            r = t/(2.0*sin(θ))
            y = sqrt(r^2-(t^2)/4.0)
            x̲ = edgeMidpoints[j].-sign(cellPressures[i]*B[i,j]+cellPressures[ii]*B[ii,j]).*(ϵ*t̂).*y
            edgeEnd1 = edgeMidpoints[j]-t̲./2.0
            triangleEdge1 = edgeEnd1.-x̲
            angle1 = atan(triangleEdge1[2],triangleEdge1[1])
            edgeEnd2 = edgeMidpoints[j]+t̲./2.0
            triangleEdge2 = edgeEnd2.-x̲
            angle2 = atan(triangleEdge2[2],triangleEdge2[1])        
            if (max(angle1, angle2) - min(angle1, angle2)) > π
                arc!(ax1, x̲, r, max(angle1, angle2), 2π+min(angle1, angle2), color=:black)
            else 
                arc!(ax1, x̲, r, angle1, angle2, color=:black)
            end
        else 
            i = findall(x->x!=0,B[:,j])[1]
            t̲ = edgeTangents[j]
            t = norm(t̲)
            t̂ = t̲./t
            k = findall(x->x<0, A[j,:])
            resultantF = 0.5*cellPressures[i]*B[i,j]*Ā[j,k].*(ϵ*edgeTangents[j])
            resultantF -= cellTensions[i]*B̄[i,j]*A[j,k].*edgeTangents[j]./edgeLengths[j]
            resultantF += 0.5*params.pressureExternal*B[i,j].*(ϵ*edgeTangents[j])
            θ = acos(normalize(resultantF)⋅t̂)
            r = t/(2.0*sin(θ))
            y = sqrt(r^2-(t^2)/4.0)
            x̲ = edgeMidpoints[j].-sign(cellPressures[i]*B[i,j]-params.pressureExternal*B[i,j]).*(ϵ*t̂).*y
            edgeEnd1 = edgeMidpoints[j]-t̲./2.0
            triangleEdge1 = edgeEnd1.-x̲
            angle1 = atan(triangleEdge1[2],triangleEdge1[1])
            edgeEnd2 = edgeMidpoints[j]+t̲./2.0
            triangleEdge2 = edgeEnd2.-x̲
            angle2 = atan(triangleEdge2[2],triangleEdge2[1])
            if (max(angle1, angle2) - min(angle1, angle2)) > π
                arc!(ax1, x̲, r, max(angle1, angle2), 2π+min(angle1, angle2), color=:black)
            else 
                arc!(ax1, x̲, r, angle1, angle2, color=:black)
            end
        end
    end
    
    # Set limits
    reset_limits!(ax1)
    recordframe!(mov)

end

save(datadir("sims", folderName, "arcs.mp4"), mov)

# pressureStrings = [@sprintf("%.2f", p) for p in cellPressures]
# annotations!(ax1, pressureStrings, Point{2,Float64}.(cellPositions))
