# Import Julia packages
using DrWatson
@quickactivate
using Revise
using LinearAlgebra
using DelimitedFiles
using SparseArrays
using StaticArrays
using CairoMakie
using UnPack
using GeometryBasics
using Random
using Colors
using JLD2

# Local modules
includet("$(projectdir())/scripts/analysisFunctions/functions.jl")

dataDirectory = "data/sims/2022-03-16-16-02-03"

function intersectionDivsCurls(dataDirectory)

    isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png")
    isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf")
    isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg")
    isdir("$dataDirectory/png") ? nothing : mkpath("$dataDirectory/png")
    isdir("$dataDirectory/pdf") ? nothing : mkpath("$dataDirectory/pdf")
    isdir("$dataDirectory/svg") ? nothing : mkpath("$dataDirectory/svg")

    # Import system data
    conditionsDict    = load("$dataDirectory/dataFinal.jld2")
    @unpack initialSystem,nVerts,nCells,nEdges,γ,λ,L₀,A₀,pressureExternal,dt,outputTotal,outputInterval,viscousTimeScale,realTimetMax,tMax,realCycleTime,nonDimCycleTime,t1Threshold = conditionsDict["params"]
    matricesDict = load("$dataDirectory/matricesFinal.jld2")
    @unpack R,A,B,Aᵀ,Ā,Āᵀ,Bᵀ,B̄,B̄ᵀ,C,cellEdgeCount,boundaryVertices,cellPositions,cellPerimeters,cellAreas,cellTensions,cellPressures,edgeLengths,edgeTangents,edgeMidpoints,F,externalF,ϵ = matricesDict["matrices"]

    # Create vector of polygons for each cell
    cellPolygons = makeCellPolygons(conditionsDict["params"],matricesDict["matrices"])

    # Find cell midpoint links T
    T = makeCellLinks(conditionsDict["params"],matricesDict["matrices"])

    # Create vector of triangles from midpoint links
    linkTriangles = makeLinkTriangles(conditionsDict["params"],matricesDict["matrices"])
    linkTriangleAreas = abs.(area.(linkTriangles))

    edgeTrapezia = makeEdgeTrapezia(conditionsDict["params"],matricesDict["matrices"])
    trapeziumAreas = abs.(area.(edgeTrapezia))

    intersections = edgeLinkMidpoints(conditionsDict["params"],matricesDict["matrices"],trapeziumAreas,T)

    q = calculateSpokes(conditionsDict["params"],matricesDict["matrices"])

    onesVec = ones(1,nCells)
    boundaryEdges = abs.(onesVec*B)
    cᵖ = boundaryEdges'.*edgeMidpoints

    vertexMidpointCurls = calculateVertexMidpointCurls(conditionsDict["params"],matricesDict["matrices"],intersections,linkTriangleAreas,q)
    vertexMidpointCurlLims = (-maximum(abs.(vertexMidpointCurls)),maximum(abs.(vertexMidpointCurls)))

    vertexMidpointDivs = calculateVertexMidpointDivs(conditionsDict["params"],matricesDict["matrices"],intersections,linkTriangleAreas,q)
    vertexMidpointDivLims = (-maximum(abs.(vertexMidpointDivs)),maximum(abs.(vertexMidpointDivs)))

    cellMidpointDivs = calculateCellMidpointDivs(conditionsDict["params"],matricesDict["matrices"],intersections,q)
    cellMidpointDivLims = (-maximum(abs.(cellMidpointDivs)),maximum(abs.(cellMidpointDivs)))

    cellMidpointCurls = calculateCellMidpointCurls(conditionsDict["params"],matricesDict["matrices"],intersections,q)
    cellMidpointCurlLims = (-maximum(abs.(cellMidpointCurls)),maximum(abs.(cellMidpointCurls)))

    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")

    # Set up figure canvas
    fig1 = Figure(resolution=(500,500))
    ax1 = Axis(fig1[1,1],aspect=DataAspect())
    hidedecorations!(ax1)
    hidespines!(ax1)

    for k=1:nVerts
        poly!(ax1,linkTriangles[k],color=[vertexMidpointCurls[k]],colorrange=vertexMidpointCurlLims,colormap=:bwr,strokecolor=(:orange,0.75),strokewidth=1)
    end
    # Plot cell polygons
    for i=1:nCells
        poly!(ax1,cellPolygons[i],color=(:white,0.0),strokewidth=1,strokecolor=(:black,1.0))
    end
    scatter!(ax1,Point2f.(intersections),color=:orange,markersize=4)
    Colorbar(fig1[1,2],limits=vertexMidpointCurlLims,colormap=:bwr)

    colsize!(fig1.layout,1,Aspect(1,1.0))
    resize_to_layout!(fig1)

    save("$dataDirectory/pdf/intersectionVertexCurls.pdf",fig1)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/intersectionVertexCurls.pdf",fig1)
    save("$dataDirectory/svg/intersectionVertexCurls.svg",fig1)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/intersectionVertexCurls.svg",fig1)
    save("$dataDirectory/png/intersectionVertexCurls.png",fig1)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/intersectionVertexCurls.png",fig1)

    # Set up figure canvas
    # fig2 = Figure(resolution=(500,500))
    # grid2 = fig2[1,1] = GridLayout()
    # ax2 = Axis(grid2[1,1][1,1],aspect=DataAspect())
    # hidedecorations!(ax2)
    # hidespines!(ax2)
    #
    # # Plot cell polygons
    # for i=1:nCells
    #     poly!(ax2,cellPolygons[i],color=(:white,0.0),strokewidth=1)
    # end
    #
    # for k=1:nVerts
    #     poly!(ax2,linkTriangles[k],color=[vertexMidpointDivs[k]],colorrange=vertexMidpointDivLims,colormap=:bwr,strokecolor=(:black,0.5),strokewidth=1)
    # end
    #
    # for j=1:nEdges
    #     edgeCells = findall(!iszero,B[:,j])
    #     if boundaryEdges[j] == 0
    #         lines!(ax2,Point2f.(cellPositions[edgeCells]),linewidth=1,color=(:white,1.0))
    #     else
    #         lines!(ax2,Point2f.([cellPositions[edgeCells[1]],cᵖ[j]]),linewidth=1,color=(:white,1.0))
    #     end
    # end
    #
    # scatter!(ax2,Point2f.(R),alpha=0.5,color=:blue)
    # scatter!(ax2,Point2f.(edgeMidpoints),alpha=0.5,color=:green)
    # scatter!(ax2,Point2f.(cellPositions),color=:red)
    # scatter!(ax2,Point2f.(intersections),color=:orange)
    # Colorbar(grid2[1,1][1,2],limits=vertexMidpointDivLims,colormap=:bwr)
    #
    # save("$dataDirectory/pdf/intersectionVertexDivs.pdf",fig2)
    # save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/intersectionVertexDivs.pdf",fig2)
    # save("$dataDirectory/svg/intersectionVertexDivs.svg",fig2)
    # save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/intersectionVertexDivs.svg",fig2)
    # save("$dataDirectory/png/intersectionVertexDivs.png",fig2)
    # save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/intersectionVertexDivs.png",fig2)

    # Set up figure canvas
    # fig3 = Figure(resolution=(500,500))
    # grid3 = fig3[1,1] = GridLayout()
    # ax3 = Axis(grid3[1,1][1,1],aspect=DataAspect())
    # hidedecorations!(ax3)
    # hidespines!(ax3)
    #
    # # Plot cell polygons
    # for i=1:nCells
    #     poly!(ax3,cellPolygons[i],color=[cellMidpointDivs[i]],colorrange=cellMidpointCurlLims,colormap=:bwr,strokewidth=1)
    # end
    #
    # for j=1:nEdges
    #     edgeCells = findall(!iszero,B[:,j])
    #     if boundaryEdges[j] == 0
    #         lines!(ax3,Point2f.(cellPositions[edgeCells]),linewidth=1,color=(:white,1.0))
    #     else
    #         lines!(ax3,Point2f.([cellPositions[edgeCells[1]],cᵖ[j]]),linewidth=1,color=(:white,1.0))
    #     end
    # end
    #
    # scatter!(ax3,Point2f.(R),alpha=0.5,color=:blue)
    # scatter!(ax3,Point2f.(edgeMidpoints),alpha=0.5,color=:green)
    # scatter!(ax3,Point2f.(cellPositions),color=:red)
    # scatter!(ax3,Point2f.(intersections),color=:orange)
    # Colorbar(grid3[1,1][1,2],limits=cellMidpointDivLims,colormap=:bwr)
    #
    # save("$dataDirectory/pdf/intersectionCellDivs.pdf",fig3)
    # save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/intersectionCellDivs.pdf",fig3)
    # save("$dataDirectory/svg/intersectionCellDivs.svg",fig3)
    # save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/intersectionCellDivs.svg",fig3)
    # save("$dataDirectory/png/intersectionCellDivs.png",fig3)
    # save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/intersectionCellDivs.png",fig3)

    # Set up figure canvas
    fig4 = Figure(resolution=(500,500))
    ax4 = Axis(fig4[1,1],aspect=DataAspect())
    hidedecorations!(ax4)
    hidespines!(ax4)

    # Plot cell polygons
    for i=1:nCells
        poly!(ax4,cellPolygons[i],color=[cellMidpointCurls[i]],colorrange=cellMidpointCurlLims,colormap=:bwr,strokewidth=1)
    end

    for k=1:nVerts
        poly!(ax4,linkTriangles[k],color=(:white,0.0),colormap=:bwr,strokecolor=(:orange,0.75),strokewidth=1)
    end

    # for j=1:nEdges
    #     edgeCells = findall(!iszero,B[:,j])
    #     if boundaryEdges[j] == 0
    #         lines!(ax4,Point2f.(cellPositions[edgeCells]),linewidth=1,color=(:white,1.0))
    #     else
    #         lines!(ax4,Point2f.([cellPositions[edgeCells[1]],cᵖ[j]]),linewidth=1,color=(:white,1.0))
    #     end
    # end

    # scatter!(ax4,Point2f.(R),alpha=0.5,color=:blue)
    # scatter!(ax4,Point2f.(edgeMidpoints),alpha=0.5,color=:green)
    # scatter!(ax4,Point2f.(cellPositions),color=:red)
    scatter!(ax4,Point2f.(intersections),color=:orange,markersize=4)
    Colorbar(fig4[1,2],limits=cellMidpointCurlLims,colormap=:bwr)
    colsize!(fig4.layout,1,Aspect(1,1.0))
    resize_to_layout!(fig4)

    save("$dataDirectory/pdf/intersectionCellCurls.pdf",fig4)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/intersectionCellCurls.pdf",fig4)
    save("$dataDirectory/svg/intersectionCellCurls.svg",fig4)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/intersectionCellCurls.svg",fig4)
    save("$dataDirectory/png/intersectionCellCurls.png",fig4)
    save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/intersectionCellCurls.png",fig4)
end
