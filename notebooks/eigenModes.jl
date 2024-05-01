# Import Julia packages
using DrWatson
using LinearAlgebra
using SparseArrays
using StaticArrays
using CairoMakie
using UnPack
using GeometryBasics
using Colors
using JLD2
using FromFile

# Local modules
@from "$(projectdir())/src/VertexModelContainers.jl" using VertexModelContainers
@from "$(projectdir())/src/AnalysisFunctions.jl" using AnalysisFunctions
@from "$(projectdir())/src/Laplacians.jl" using Laplacians

folderName = "L₀=0.75_nCells=61_pressureExternal=0.5_realTimetMax=86400.0_stiffnessFactor=1.0_γ=0.1_24-05-01-16-02-18"

files = [datadir("sims", folderName, "frameData", f) for f in readdir(datadir("sims", folderName, "frameData")) if occursin(".jld2",f)]
@unpack R, matrices, params = load(files[end]; 
        typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer"=>ParametersContainer, 
        "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer"=>MatricesContainer))
@unpack nVerts, nCells = params

T = makeCellLinks(params,matrices)

edgeTrapezia = makeEdgeTrapezia(R,params,matrices)
trapeziumAreas = abs.(area.(edgeTrapezia))

linkTriangles = makeLinkTriangles(R, params,matrices)

cellPolygons = makeCellPolygons(R, params,matrices)

Lf = makeLf(params,matrices,trapeziumAreas)
decompositionLf = (eigen(Matrix(Lf))).vectors

Lᵥ = makeLv(params,matrices,abs.(area.(linkTriangles)),trapeziumAreas)
decompositionLv = (eigen(Matrix(Lᵥ))).vectors

# Set up figure canvas
fig = Figure(size=(1500,1600))

# for x=1:5
#     for y=1:4
#         println("$(((y-1)*5 + x)+1), Lv")
#         eigenvectorIndex = ((y-1)*5 + x)+1
#         lims = (-maximum(abs.(decompositionLv[:,eigenvectorIndex])),maximum(abs.(decompositionLv[:,eigenvectorIndex])))
#         ax = Axis(grid1[y,x],aspect=DataAspect())
#         hidedecorations!(ax)
#         hidespines!(ax)
#         for k=1:nVerts
#             poly!(ax,linkTriangles[k],color=[decompositionLv[k,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
#         end
#         for i=1:nCells
#             poly!(ax,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
#         end
#         Label(grid1[y,x,Bottom()],
#                 L"k=%$eigenvectorIndex",
#                 fontsize = 16,
#         )
#     end
# end
for x=1:5
    for y=1:4
        println("$(((y-1)*5 + x)+(nVerts-20)), Lv")
        eigenvectorIndex = ((y-1)*5 + x)+(nVerts-20)
        lims = (-maximum(abs.(decompositionLv[:,eigenvectorIndex])),maximum(abs.(decompositionLv[:,eigenvectorIndex])))
        # ax = Axis(fig[y+4,x],aspect=DataAspect())
        ax = Axis(fig[y,x],aspect=DataAspect())
        hidedecorations!(ax)
        hidespines!(ax)
        for k=1:nVerts
            poly!(ax,linkTriangles[k],color=decompositionLv[k,eigenvectorIndex],colorrange=lims,colormap=:bwr,strokewidth=0,strokecolor=(:black,0.25)) #:bwr
        end
        for i=1:nCells
            if matrices.μ[i] < 1.5
                poly!(ax,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,0.1),strokewidth=0.1)
            end
        end
        for i=1:nCells
            if matrices.μ[i] > 1.5
                poly!(ax,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=0.1)
            end
        end
        Label(grid1[y+4,x,Bottom()],
                L"k=%$eigenvectorIndex",
                fontsize = 16,
        )
    end
end
# Label(grid1[9,:],
#         L"(a)",
#         fontsize = 32,
# )


# vlineAx = Axis(superGrid[1,2])
# hidedecorations!(vlineAx)
# hidespines!(vlineAx)
# vlines!(vlineAx,0.0,color=:black)
# colsize!(superGrid, 2, Aspect(1, 0.01))


# signInversions = [3, 5, 9, 10, 11, 12, 13, 17, 19, 20, 21]

# grid2 = superGrid[1,3] = GridLayout()
# for x=1:5
#     @show x
#     for y=1:4
#         @show y
#         eigenvectorIndex = ((y-1)*5 + x)+1
#         lims = (-maximum(abs.(decompositionLf[:,eigenvectorIndex])),maximum(abs.(decompositionLf[:,eigenvectorIndex])))
#         ax = Axis(grid2[y,x],aspect=DataAspect())
#         hidedecorations!(ax)
#         hidespines!(ax)
#         # if eigenvectorIndex in signInversions
#         #     for i=1:nCells
#         #         poly!(ax,cellPolygons[i],color=[-decompositionLf[i,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokecolor=(:black,1.0),strokewidth=1) #:bwr
#         #     end
#         # else
#             for i=1:nCells
#                 poly!(ax,cellPolygons[i],color=[decompositionLf[i,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokecolor=(:black,1.0),strokewidth=1) #:bwr
#             end
#         # end
#         Label(grid2[y,x,Bottom()],
#                 L"i=%$eigenvectorIndex",
#                 fontsize = 16,
#         )
#     end
# end
for x=1:5
    for y=1:4
        println("$(((y-1)*5 + x)+(nVerts-20)), Lf")
        eigenvectorIndex = ((y-1)*5 + x)+(nCells-20)
        lims = (-maximum(abs.(decompositionLf[:,eigenvectorIndex])),maximum(abs.(decompositionLf[:,eigenvectorIndex])))
        ax = Axis(fig[y+4,x],aspect=DataAspect())
        hidedecorations!(ax)
        hidespines!(ax)
        for i=1:nCells
            if matrices.μ[i] < 1.5
                poly!(ax,cellPolygons[i],color=decompositionLf[i,eigenvectorIndex],colorrange=lims,colormap=:bwr,strokecolor=(:black,0.1),strokewidth=0.1)
            end
        end
        for i=1:nCells
            if matrices.μ[i] > 1.5
                poly!(ax,cellPolygons[i],color=decompositionLf[i,eigenvectorIndex],colorrange=lims,colormap=:bwr,strokecolor=(:black,1.0),strokewidth=0.1)
            end
        end
        # Label(grid2[y+4,x,Bottom()],
        #         L"i=%$eigenvectorIndex",
        #         fontsize = 16,
        # )
    end
end
# Label(grid2[9,:],
#         L"(b)",
#         fontsize = 32,
# )

# resize_to_layout!(fig)

display(fig)
save(datadir("sims", folderName,"eigenModes.png"),fig)
