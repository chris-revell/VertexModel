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

@from "TopologyChange.jl" using TopologyChange
@from "Model.jl" using Model
@from "EdgeAblation" using EdgeAblation

# Pick the equilibrated system: 
# dateString = ""
# parameterSetLabel = "" 

# !isdir(datadir("multipleRuns", dateString, parameterSetLabel, "ablationLoop")) ? mkpath(datadir("multipleRuns", dateString, paramterSetLabel, "ablationLoop")) : nothing 

# dataDict = load(datadir("multipleRuns",dateString, parameterSetLabel, "$(parameterSetLabel)_equilibriumPhase.jld2");
#                     typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer" => MatricesContainer, 
#                                 "VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer" => ParametersContainer
#                     )
#                 )

dataDict = load(datadir("/Users/charlietaylorbarca/Documents/GitHub/VertexModel/data/sims/charlie-free-boundaries/(VIII)/Growth and equilibrium/26-06-09-14-03-05_nCells=795_Λ_AA=-0.1_Λ_AB=-0.2_Λ_BB=-0.3_β=0.0_γ=0.05/frameData/systemData099.jld2");
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
        energyModel = params 
@unpack Ā,
        B̄,
        cellTensions,
        cellPressures,
        edgeLengths,
        edgeTangents,
        F,
        externalF,
        ϵ,
        boundaryVertices,
        boundaryEdges,
        vertexAreas,
        Λs = matrices

for jAblated = 1

    if jAblated in boundaryEdges
        break
    end

    # Find the vertices at either end of the edge: 
    params.k_tracked = findall(x -> x!=0, @view matrices.A[params.jAblated,:])

    edgeAblation!(jAblated, params, matrices)
    topologyChange!(R,params,matrices)
    spatialData!(R, params, matrices)

    # Calculate the resultant force at each k_tracked after ablation: 
    F = zeros(SVector{2,Float64}, 2, nCells)
    dR = zeros(SVector{2,Float64}, 2)
    for k in k_tracked
        for j in nzrange(A, k) # iterate over the nonzero entries for vertex k 
            for i in nzrange(B, rowvals(A)[j]) # rowvals(A) gives the row indices of nonzero entries of A
                
                # Force components from cell pressure perpendicular to edge tangents - the area derivative wrt. vertex position of Energy from pressure
                F[k, rowvals(B)[i]] += 0.5 * cellPressures[rowvals(B)[i]] * B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ * edgeTangents[rowvals(A)[j]])
                # Force components from cell membrane tension parallel to edge tangents 
                F[k, rowvals(B)[i]] -= cellTensions[rowvals(B)[i]] * B̄[rowvals(B)[i], rowvals(A)[j]] * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
                # Force on vertex from external pressure -- only applies to boundary vertices 
                if boundaryType == "free"
                    externalF[k] += boundaryVertices[k] * (0.5 * pressureExternal * B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ * edgeTangents[rowvals(A)[j]])) # 0 unless boundaryVertices != 0
                end
                
                if energyModel == "quadratic2pops"
                    # We have the separate edge tension term in that case: 
                    # Factor of 1/2 because this is added for each cell that meets at j
                    F[k, rowvals(B)[i]] -=  0.5 * Λs[rowvals(A)[j]] * B̄[rowvals(B)[i], rowvals(A)[j]] *  A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
                end 
                
            end
            # Force on vertex from peripheral tension -- only for boundary edges 
            if boundaryType == "free"
                externalF[k] -= boundaryEdges[rowvals(A)[j]] * peripheralTension * (peripheryLength - sqrt(π * nCells)) * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
            end
        end

        dR[k] = sum(F[k, :])

    end
    println(dR)



end

jldsave(datadir("multipleRuns", dateString, parameterSetLabel, "ablationLoop"); )

