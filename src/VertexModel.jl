#
#  VertexModel.jl
#  VertexModel
#
#  Created by Christopher Revell on 31/01/2021.
#
#

module VertexModel

# Julia packages
using LinearAlgebra
using DelimitedFiles
using SparseArrays
using StaticArrays
using LoopVectorization
using Plots
using UnPack

# Local modules
include("TopologyChange.jl"); using .TopologyChange
include("CreateRunDirectory.jl"); using .CreateRunDirectory
include("Visualise.jl"); using .Visualise
include("Initialise.jl"); using .Initialise
include("Iterate.jl"); using .Iterate


function vertexModel(initialSystem,realTimetMax,γ,λ,tStar,dt,preferredArea,pressureExternal,outputTotal,t1Threshold,outputToggle)

    #BLAS.set_num_threads(4)

    # Input parameters
    # realTimetMax     (eg. 10000.0)  Real time maximum system run time /seconds
    # γ                (eg. 0.2    )  Parameters in energy relaxation
    # λ                (eg. -0.3   )  Parameters in energy relaxation
    # tStar            (eg. 20.0   )  Relaxation rate, approx from Sarah's data.
    # dt               (eg. 0.01   )  Non dimensionalised time step
    # preferredArea    (eg. 1.0    )  Cell preferred area (1.0 by default)
    # pressureExternal (eg. 0.1    )  External pressure applied isotropically to system boundary
    # outputTotal      (eg. 20     )  Number of data outputs
    # t1Threshold      (eg. 0.01   )  Edge length at which a T1 transition is triggered

    R,tempR,ΔR,nVerts,nCells,nEdges,γ,λ,preferredPerimeter,preferredArea,pressureExternal,dt,outputInterval,tStar,realTimetMax,tMax,t1Threshold,A,B,Aᵀ,Ā,Āᵀ,Bᵀ,B̄,B̄ᵀ,C,cellEdgeCount,boundaryVertices,cellPositions,cellPerimeters,cellOrientedAreas,cellAreas,cellTensions,cellPressures,edgeLengths,edgeTangents,edgeMidpoints,vertexEdges,vertexCells,F,ϵ = initialise(initialSystem,realTimetMax,γ,λ,preferredArea,pressureExternal,dt,tStar,outputTotal,t1Threshold)

    # Create output directory in which to store results and parameters, and initialise animation object for plotting
    if outputToggle==1
        folderName = createRunDirectory(R,nVerts, nCells, nEdges, γ, λ, preferredPerimeter, preferredArea, pressureExternal, dt, outputInterval, tStar, realTimetMax, tMax,A,B)
        anim = Animation()
    end

    # Initial setup of matrices based on system topology
    topologyChange!(A,B,Aᵀ,Ā,Āᵀ,Bᵀ,B̄,B̄ᵀ,C,cellEdgeCount,boundaryVertices)

    t = 1E-8   # Initial time is very small but slightly above 0 to avoid issues with remainders in output interval calculation

    while t<tMax

        # 4 step Runge-Kutta integration
        # 1st step of Runge-Kutta
        iterate1!(6.0,R,ΔR,F,dt,nVerts,nCells,nEdges,γ,preferredPerimeter,preferredArea,A,B,Ā,B̄,C,cellEdgeCount,cellPositions,cellPerimeters,cellOrientedAreas,cellAreas,cellTensions,cellPressures,edgeLengths,edgeTangents,edgeMidpoints,boundaryVertices,ϵ,t1Threshold)

        # 2nd step of Runge-Kutta
        iterate2!(2.0,3.0,R,tempR,ΔR,F,dt,nVerts,nCells,nEdges,γ,preferredPerimeter,preferredArea,A,B,Ā,B̄,C,cellEdgeCount,cellPositions,cellPerimeters,cellOrientedAreas,cellAreas,cellTensions,cellPressures,edgeLengths,edgeTangents,edgeMidpoints,boundaryVertices,ϵ)

        # 3rd step of Runge-Kutta
        iterate2!(2.0,3.0,R,tempR,ΔR,F,dt,nVerts,nCells,nEdges,γ,preferredPerimeter,preferredArea,A,B,Ā,B̄,C,cellEdgeCount,cellPositions,cellPerimeters,cellOrientedAreas,cellAreas,cellTensions,cellPressures,edgeLengths,edgeTangents,edgeMidpoints,boundaryVertices,ϵ)

        # 4th step of Runge-Kutta
        iterate2!(1.0,6.0,R,tempR,ΔR,F,dt,nVerts,nCells,nEdges,γ,preferredPerimeter,preferredArea,A,B,Ā,B̄,C,cellEdgeCount,cellPositions,cellPerimeters,cellOrientedAreas,cellAreas,cellTensions,cellPressures,edgeLengths,edgeTangents,edgeMidpoints,boundaryVertices,ϵ)


        # Result of Runge-Kutta steps
        R .+= ΔR
        t +=dt

        if t%outputInterval<dt && outputToggle==1
            visualise(anim,R,nEdges,nVerts,nCells,A,Ā,B̄,C,F,cellPositions,edgeTangents,edgeMidpoints,boundaryVertices,vertexEdges)
            println("$t/$tMax)")
        end

    end

    # If outputToggle==1, save animation object as an animated gif
    outputToggle==1 ? gif(anim, "data/sims/$folderName/animated.gif", fps = 10) : nothing

end

export vertexModel

end
