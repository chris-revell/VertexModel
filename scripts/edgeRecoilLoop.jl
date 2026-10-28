# Script to run over edges of a system and calculate the speed of each vertex 

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

include(srcdir("TopologyChange.jl")); using .TopologyChange
include(srcdir("Model.jl")); using .Model
include(srcdir("EdgeAblation.jl")); using .EdgeAblation
include(srcdir("SpatialData.jl")); using .SpatialData
include(srcdir("AnalysisFunctions.jl")); using .AnalysisFunctions

# Pick the equilibrated system: 
dateString = "26-09-10-10-49-26"



# dataDict = load(datadir("multipleRuns",dateString, parameterSetLabel, "$(parameterSetLabel)_equilibriumPhase.jld2");
#                     typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer" => MatricesContainer, 
#                                 "VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer" => ParametersContainer
#                     )
#                 )

# dataDict = load("C:/Users/z28439ct/Documents/GitHub/VertexModel/data/multipleRuns/26-09-10-10-49-26/(I)/26-09-12-14-42-03_nCells=1273_Λ_AA=-0.2_Λ_AB=-0.2_Λ_BB=-0.2_β=0.0_γ=0.05/frameData/systemData099.jld2";
#                 typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer" => MatricesContainer, 
#                             "VertexModel...VertexModelContainers.jl.VertexModelContainers.ParametersContainer" => ParametersContainer))


jld2pathVec = ["data/multipleRuns/26-09-10-10-49-26/(O)/26-09-12-14-42-03_nCells=1273_Λ_AA=-0.2_Λ_AB=-0.2_Λ_BB=-0.2_β=0.0_γ=0.05/frameData/systemData099.jld2", #(O)
                "data/multipleRuns/26-09-10-10-49-26/(P)/26-09-12-13-53-32_nCells=1273_Λ_AA=-0.1_Λ_AB=-0.2_Λ_BB=-0.25_β=0.0_γ=0.05/frameData/systemData099.jld2",#(P)
                "data/multipleRuns/26-09-10-10-49-26/(Q)/26-09-12-18-25-13_nCells=1273_Λ_AA=-0.1_Λ_AB=-0.2_Λ_BB=-0.3_β=0.0_γ=0.05/frameData/systemData099.jld2",#(Q)
                "data/multipleRuns/26-09-10-10-49-26/(R)/26-09-12-13-39-48_nCells=1273_Λ_AA=-0.15_Λ_AB=-0.2_Λ_BB=-0.3_β=0.0_γ=0.05/frameData/systemData099.jld2",#(R)
                "data/multipleRuns/26-09-10-10-49-26/(S)/26-09-13-03-46-48_nCells=1169_Λ_AA=-0.3_Λ_AB=-0.2_Λ_BB=-0.3_β=0.0_γ=0.05/frameData/systemData099.jld2",#(S)
                "data/multipleRuns/26-09-10-10-49-26/(T)/26-09-12-14-02-50_nCells=1273_Λ_AA=-0.3_Λ_AB=-0.2_Λ_BB=-0.15_β=0.0_γ=0.05/frameData/systemData099.jld2",#(T)
                "data/multipleRuns/26-09-10-10-49-26/(U)/26-09-12-14-08-47_nCells=1273_Λ_AA=-0.3_Λ_AB=-0.2_Λ_BB=-0.1_β=0.0_γ=0.05/frameData/systemData099.jld2",#(U)
                "data/multipleRuns/26-09-10-10-49-26/(V)/26-09-12-13-57-33_nCells=1273_Λ_AA=-0.25_Λ_AB=-0.2_Λ_BB=-0.1_β=0.0_γ=0.05/frameData/systemData099.jld2",#(V)
                "data/multipleRuns/26-09-10-10-49-26/(KRAS V12 TEST 1)/26-10-21-09-55-13_nCells=1273_Λ_AA=-0.325_Λ_AB=-0.2_Λ_BB=-0.075_β=0.0_γ=0.05/frameData/systemData072.jld2",]#((KRAS V12 TEST 1))


parameterLabelVec = ["(O)","(P)","(Q)","(R)","(S)","(T)","(U)","(V)","(KRAS V12 TEST 1)"]

# Flag for cases where the recoil data has been calculated already: 
importRecoilData = false

for (simInd,jld2path) in enumerate(jld2pathVec)


    dataDict = load(jld2path;
                    typemap=Dict("VertexModel...\\VertexModelContainers.jl.VertexModelContainers.MatricesContainer" => MatricesContainer, 
                                "VertexModel...\\VertexModelContainers.jl.VertexModelContainers.ParametersContainer" => ParametersContainer)
                                )

    parameterSetLabel = parameterLabelVec[simInd]

    !isdir(datadir("multipleRuns", dateString, parameterSetLabel, "ablationLoop")) ? mkpath(datadir("multipleRuns", dateString, parameterSetLabel, "ablationLoop")) : nothing 


    # Import system data
    @unpack R, params, matrices = dataDict
    @unpack nVerts,
            nCells,
            nEdges,
            pressureExternal,
            peripheralTension,
            vertexWeighting,
            energyModel,
            boundaryType,
            k_tracked = params 
    @unpack A,
            B,
            Ā,
            B̄,
            cellTensions,
            cellPressures,
            edgeLengths,
            edgeTangents,
            ϵ,
            boundaryVertices,
            boundaryEdges,
            vertexAreas,
            Λs,
            jsAfterAblation,
            edgeLabels = matrices

    # Store deep copies of the incidence and vertex matrices since ablation function will modify them in place: 
    R_orig = deepcopy(R)
    params_orig = deepcopy(params)
    matrices_orig = deepcopy(matrices)

    if !importRecoilData

        # Initialise matrices to store speed and rotation on each edge if it is ablated 
        edgeSpeeds = zeros(Float64, nEdges)
        edgeRotations = zeros(Float64, nEdges)

        for jAblated = 1:nEdges
            println("j = ", jAblated, "/$nEdges")
        # jAblated = 1

            if matrices_orig.boundaryEdges[jAblated] == 1
                continue
            end

            #

            R_l = deepcopy(R_orig)
            params_l = deepcopy(params_orig)
            matrices_l = deepcopy(matrices_orig)

            # Find the vertices at either end of the edge: 
            params_l.k_tracked = findall(x -> x!=0, @view matrices_l.A[jAblated,:])


            EdgeAblation.edgeAblation!(jAblated, params_l, matrices_l)
            TopologyChange.topologyChange!(R_l, params_l, matrices_l)
            SpatialData.spatialData!(R_l, params_l, matrices_l)

            # Unpack deep copies of matrices: 
            A_l, B_l, Ā_l, B̄_l = matrices_l.A, matrices_l.B, matrices_l.Ā, matrices_l.B̄
            cellTensions_l, cellPressures_l = matrices_l.cellTensions, matrices_l.cellPressures
            edgeLengths_l, edgeTangents_l   = matrices_l.edgeLengths, matrices_l.edgeTangents
            ϵ_l, Λs_l = matrices_l.ϵ, matrices_l.Λs
            energyModel,k_tracked_l = params_l.energyModel,params_l.k_tracked

            # Calculate the resultant force at each k_tracked after ablation: 
            ablatedF = zeros(SVector{2,Float64}, 2, nCells)
            dR = zeros(SVector{2,Float64}, 2)

            for kInd in enumerate(k_tracked_l)
                k = kInd[2] # vertex index
                kInd = kInd[1] # index in k_tracked array
                for j in nzrange(A_l, k) # iterate over the nonzero entries for vertex k 
                    for i in nzrange(B_l, rowvals(A_l)[j]) # rowvals(A) gives the row indices of nonzero entries of A
                        
                        # Force components from cell pressure perpendicular to edge tangents - the area derivative wrt. vertex position of Energy from pressure
                        ablatedF[kInd, rowvals(B_l)[i]] += 0.5 * cellPressures_l[rowvals(B_l)[i]] * B_l[rowvals(B_l)[i], rowvals(A_l)[j]] * Ā_l[rowvals(A_l)[j], k] .* (ϵ * edgeTangents_l[rowvals(A_l)[j]])
                        # Force components from cell membrane tension parallel to edge tangents 
                        ablatedF[kInd, rowvals(B_l)[i]] -= cellTensions_l[rowvals(B_l)[i]] * B̄_l[rowvals(B_l)[i], rowvals(A_l)[j]] * A_l[rowvals(A_l)[j], k] .* edgeTangents_l[rowvals(A_l)[j]] ./ edgeLengths_l[rowvals(A_l)[j]]
                        # Force on vertex from external pressure -- only applies to boundary vertices 
                        # if boundaryType == "free"
                        #     externalF[kInd] += boundaryVertices[k] * (0.5 * pressureExternal * B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ * edgeTangents[rowvals(A)[j]])) # 0 unless boundaryVertices != 0
                        # end
                        
                        if energyModel == "quadratic2pops"
                            # We have the separate edge tension term in that case: 
                            # Factor of 1/2 because this is added for each cell that meets at j
                            ablatedF[kInd, rowvals(B_l)[i]] -=  0.5 * Λs_l[rowvals(A_l)[j]] * B̄_l[rowvals(B_l)[i], rowvals(A_l)[j]] *  A_l[rowvals(A_l)[j], k] .* edgeTangents_l[rowvals(A_l)[j]] ./ edgeLengths_l[rowvals(A_l)[j]]
                        end 
                        
                    end
                    # Force on vertex from peripheral tension -- only for boundary edges 
                    # if boundaryType == "free"
                    #     externalF[kInd] -= boundaryEdges[rowvals(A)[j]] * peripheralTension * (peripheryLength - sqrt(π * nCells)) * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
                    # end
                end

                
                dR[kInd] = sum(ablatedF[kInd, :])
                # Scale with vertex area - drag term
                dR[kInd] ./ vertexAreas[k]

            end

            # println("dR = ",dR)
            dR̄= sum(dR)/2
            dR̂ = dR[1] - dR̄
            Δr = R_orig[k_tracked_l[1]] - R_orig[k_tracked_l[2]]

            edgeSpeeds[jAblated] = norm(dR̂)
            edgeRotations[jAblated] = Δr[1]*dR̂[2] - Δr[2]*dR̂[1] # take the cross produce with the edge. Positive = anticlockwise; negative = clockwise 

            # println("edge $jAblated: speed=$(edgeSpeeds[jAblated]), rotation=$(edgeRotations[jAblated])")

        end

    else

        # Recoil data has already been computed - extract rotations and speeds
        df = load(datadir("multipleRuns", dateString, parameterSetLabel,"ablationLoop","ablationData.jld2"))
        
        edgeSpeeds = df["edgeSpeeds"]
        edgeRotations = df["edgeRotations"]
    end

    # jldsave(datadir("multipleRuns", dateString, parameterSetLabel, "ablationLoop.jld2"); dR)
    # jldsave()

    # Now plot the speed and rotation of each edge: 

    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    edgeAblationFig = Figure(size=(1200,1200))
    grid = edgeAblationFig[1,1] = GridLayout()
    rotationAx = Axis(grid[1,1],aspect=DataAspect())
    speedAx = Axis(grid[1,2],aspect=DataAspect())


    # rotationAx.title = "Edge rotation after ablation"
    # speedAx.title = "Edge speed after ablation"
        

    hidedecorations!(rotationAx)
    hidedecorations!(speedAx)
    hidespines!(rotationAx)
    hidespines!(speedAx)

    interfaceBoundaryEdges = findall(x -> x==2,matrices_orig.edgeLabels)

    # Calculate edge vectors: 
    edgeVectors = Point2f[]
    for j in 1:nEdges
        verts = findall(x -> x != 0, @view matrices_orig.A[j, :])
        push!(edgeVectors, Point2f(R_orig[verts[1]]), Point2f(R_orig[verts[2]]))
    end

    xs = [p[1] for p in edgeVectors]
    ys = [p[2] for p in edgeVectors]
    xlims!(rotationAx, minimum(xs), maximum(xs))
    ylims!(rotationAx, minimum(ys), maximum(ys))
    xlims!(speedAx, minimum(xs), maximum(xs))
    ylims!(speedAx, minimum(ys), maximum(ys))

    cmapRotation = cgrad([
            RGB(0.0, 0.0, 1.0),    # blue
            RGB(1.0, 1.0, 1.0),   # white, zero
            RGB(1.0, 0.0, 0.0)   # red
        ], 256)
    climsRotation = (-maximum(abs.(edgeRotations)), maximum(abs.(edgeRotations)))
    # Colour bar exclusing exterior vertices 

    # cmapSpeed = cgrad([
    #         RGB(1.0, 1.0, 1.0),   # 0%   White
    #         RGB(1.0, 0.8, 0.9),   # 25%  Light Pink
    #         RGB(1.0, 0.0, 0.5),   # 50%  Hot Pink
    #         RGB(0.6, 0.0, 0.8),   # 75%  Purple
    #         RGB(0.3, 0.0, 0.5)    # 100% Dark Purple
    #     ])
    cmapSpeed = cgrad([
                RGB(1,1,1),       # white
                RGB(1,1,1),       # keep white longer
                RGB(1,0,1),       # magenta
                RGB(0.4,0.1,0.5), # purple
                RGB(0,0.85,0.85),  # teal
                RGB(0,0.85,0.85)  # teal
            ])
    climsSpeed = (0.0, maximum(edgeSpeeds))

    linesegments!(rotationAx, edgeVectors; color = edgeRotations,colorrange = climsRotation, colormap = cmapRotation, linewidth=3)
    linesegments!(speedAx, edgeVectors; color = edgeSpeeds,colorrange = climsSpeed, colormap = cmapSpeed, linewidth=3)

    # cbarRotation = Colorbar(grid[2,1],colormap = cmapRotation, colorrange=climsRotation, label="Edge rotation", width=20,height=Relative(0.6))
    # cbarSpeed = Colorbar(grid[2,2],colormap = cmapSpeed, colorrange=climsSpeed, label="Edge Speed", width=20,height=Relative(0.6))

    # # rowsize!(grid, 2, Fixed(60))   # colorbars get a fixed 60px strip; row 1 takes the rest automatically
    # rowsize!(grid, 1, Relative(0.9))
    # rowsize!(grid, 2, Relative(0.1))

    save(datadir("multipleRuns", dateString, parameterSetLabel,"ablationLoop", "ablationFigure.png"),edgeAblationFig)
    # save(datadir("multipleRuns", dateString, parameterSetLabel,"ablationLoop", "colourBar.png"),edgeAblationFig)

    jldsave(datadir("multipleRuns", dateString, parameterSetLabel,"ablationLoop","ablationData.jld2");edgeSpeeds,edgeRotations)

end