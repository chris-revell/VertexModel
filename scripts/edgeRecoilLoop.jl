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

# Pick the equilibrated system: 
dateString = "26-09-11-16-18-04"
parameterSetLabel = "(I)" 

!isdir(datadir("multipleRuns", dateString, parameterSetLabel, "ablationLoop")) ? mkpath(datadir("multipleRuns", dateString, parameterSetLabel, "ablationLoop")) : nothing 

# dataDict = load(datadir("multipleRuns",dateString, parameterSetLabel, "$(parameterSetLabel)_equilibriumPhase.jld2");
#                     typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer" => MatricesContainer, 
#                                 "VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer" => ParametersContainer
#                     )
#                 )

dataDict = load("/Users/user/The University of Manchester Dropbox/Charlotte Taylor Barca/JULIA/VertexModel/data/sims/charlie-free-boundaries/26-09-12-10-18-32_nCells=91_Λ_AA=-0.2_Λ_AB=-0.2_Λ_BB=-0.2_β=0.0_γ=0.05/frameData/systemData037.jld2";
                typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer" => MatricesContainer, 
                            "VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer" => ParametersContainer))

# Import system data
@unpack R, params, matrices = dataDict
@unpack nVerts,
        nCells,
        nEdges,
        pressureExternal,
        peripheralTension,
        vertexWeighting,
        energyModel,
        boundaryType = params 
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
        jsAfterAblation = matrices

# for jAblated = 1
jAblated = 1

    # if jAblated in boundaryEdges
    #     break
    # end

    # Find the vertices at either end of the edge: 
    params.k_tracked = findall(x -> x!=0, @view matrices.A[jAblated,:])

    EdgeAblation.edgeAblation!(jAblated, params, matrices)
    TopologyChange.topologyChange!(R,params,matrices)
    SpatialData.spatialData!(R, params, matrices)

    # Calculate the resultant force at each k_tracked after ablation: 
    ablatedF = zeros(SVector{2,Float64}, 2, nCells)
    dR = zeros(SVector{2,Float64}, 2)
    for kInd in enumerate(params.k_tracked)
        k = kInd[2] # vertex index
        kInd = kInd[1] # index in k_tracked array
        for j in nzrange(A, k) # iterate over the nonzero entries for vertex k 
            for i in nzrange(B, rowvals(A)[j]) # rowvals(A) gives the row indices of nonzero entries of A
                
                # Force components from cell pressure perpendicular to edge tangents - the area derivative wrt. vertex position of Energy from pressure
                ablatedF[kInd, rowvals(B)[i]] += 0.5 * cellPressures[rowvals(B)[i]] * B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ * edgeTangents[rowvals(A)[j]])
                # Force components from cell membrane tension parallel to edge tangents 
                ablatedF[kInd, rowvals(B)[i]] -= cellTensions[rowvals(B)[i]] * B̄[rowvals(B)[i], rowvals(A)[j]] * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
                # Force on vertex from external pressure -- only applies to boundary vertices 
                # if boundaryType == "free"
                #     externalF[kInd] += boundaryVertices[k] * (0.5 * pressureExternal * B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ * edgeTangents[rowvals(A)[j]])) # 0 unless boundaryVertices != 0
                # end
                
                if energyModel == "quadratic2pops"
                    # We have the separate edge tension term in that case: 
                    # Factor of 1/2 because this is added for each cell that meets at j
                    ablatedF[kInd, rowvals(B)[i]] -=  0.5 * Λs[rowvals(A)[j]] * B̄[rowvals(B)[i], rowvals(A)[j]] *  A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
                end 
                
            end
            # Force on vertex from peripheral tension -- only for boundary edges 
            # if boundaryType == "free"
            #     externalF[kInd] -= boundaryEdges[rowvals(A)[j]] * peripheralTension * (peripheryLength - sqrt(π * nCells)) * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
            # end
        end

        dR[kInd] = sum(ablatedF[kInd, :])

        

    end
    println("dR = ",dR)
    dR̄= sum(dR)/2
    dR̂ = dR[1] - dR̄

    println("dR̄ = ",dR̄)
    println("dR̂ = ",dR̂)



# end

# jldsave(datadir("multipleRuns", dateString, parameterSetLabel, "ablationLoop.jld2"); dR)
# jldsave()
