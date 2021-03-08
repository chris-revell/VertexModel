#
#  Simulate.jl
#  VertexModelJL
#
#  Created by Christopher Revell on 31/01/2021.
#
#

module Simulate

# Julia packages
using LinearAlgebra
using DelimitedFiles
using SparseArrays

# Local modules
using TopologyChange
using CreateRunDirectory
using SpatialData
using CalculateForce
using SingleHexagon
using Visualise

@inline @views function simulate(initialSystem)

    # Parameters
    realTimetMax = 200.0     # Real time maximum system run time /seconds
    gamma        = 0.172     # Parameters in energy relaxation. Hard wired from data.
    lamda        = -0.259    # Parameters in energy relaxation. Hard wired from data.
    tStar        = 20.0      # Relaxation rate. Approx from Sarah's data.
    dt           = 0.01      # Non dimensionalised time step
    ϵ            = [0.0 1.0
                   -1.0 0.0] # Antisymmetric rotation matrix

    # Derived parameters
    # NB Preferred area = 1.0 by default
    tMax               = realTimetMax/tStar  # Non dimensionalised maximum system run time
    outputInterval     = tMax/10.0           # Time interval for storing system data (non dimensionalised)
    preferredPerimeter = -lamda/(2*gamma)    # Cell preferred perimeter
    preferredArea      = 1.0                 # Cell preferred area

    if initialSystem=="single"
        A,B,R = singleHexagon()
    else
        # Import system matrices from file
        A = sparse(readdlm("input/$(initialSystem)_A.txt",',',Int64,'\n')) # Incidence matrix. Rows => edges; columns => vertices.
        B = sparse(readdlm("input/$(initialSystem)_B.txt",',',Int64,'\n')) # Incidence matrix. Rows => cells; columns => edges. Values +/-1 for orientation
        R = readdlm("input/$(initialSystem)_R.txt",',',Float64,'\n')       # Coordinates of vertices
    end

    # Infer system information from matrices
    nCells            = size(B)[1]                   # Number of cells
    nEdges            = size(A)[1]                   # Number of edges
    nVerts            = size(A)[2]                   # Number of vertices
    # Preallocate system arrays
    Aᵀ                = spzeros(Int64,nVerts,nEdges) # Transpose of incidence matrix A
    Ā                 = spzeros(Int64,nEdges,nVerts) # Undirected adjacency matrix from absolute values of incidence matrix A
    Āᵀ                = spzeros(Int64,nVerts,nEdges) # Undirected adjacency matrix from absolute values of transpose of incidence matrix Aᵀ
    Bᵀ                = spzeros(Int64,nEdges,nCells) # Transpose of incidence matrix B
    B̄                 = spzeros(Int64,nCells,nEdges) # Undirected adjacency matrix from absolute values of incidence matrix B
    B̄ᵀ                = spzeros(Int64,nEdges,nCells) # Undirected adjacency matrix from absolute values of transpose of incidence matrix Bᵀ
    C                 = spzeros(Int64,nCells,nVerts) # C adjacency matrix. Rows => cells; Columns => vertices. = 0.5*B̄*Ā
    tempR             = zeros(nVerts,2)              # Array to store temporary positions in Runge-Kutta integration
    ΔR                = zeros(nVerts,2)              # Array to store change in R during Runge-Kutta integration
    cellEdgeCount     = zeros(Int64,nCells,1)        # 1D matrix containing number of edges around each cell, found by summing columns of B̄
    boundaryVertices  = zeros(Int64,nVerts,1)        # 1D matrix containing labels of vertices at system boundary
    cellPositions     = zeros(nCells,2)              # 2D matrix of cell centre positions
    cellPerimeters    = zeros(nCells,1)              # 1D matrix of scalar cell perimeter lengths
    cellOrientedAreas = zeros(nCells,2,2)            # 3D array of oriented cell areas. Each row is a 2x2 antisymmetric matrix of the form [0 A / -A 0] where A is the scalar cell area
    cellAreas         = zeros(nCells,1)              # 1D matrix of scalar cell areas
    cellTensions      = zeros(nCells,1)              # 1D matrix of scalar cell tensions
    cellPressures     = zeros(nCells,1)              # 1D matrix of scalar cell pressures
    edgeLengths       = zeros(nEdges,1)              # 1D matrix of scalar edge lengths
    edgeTangents      = zeros(nEdges,2)              # 2D matrix of tangent vectors for each edge (magnitude = edge length)
    edgeMidpoints     = zeros(nEdges,2)              # 2D matrix of position coordinates for each edge midpoint
    F                 = zeros(nVerts,nCells,2)       # 3D array containing force vectors on vertex k from cell i, Fᵢₖ

    # Create output directory in which to store results and parameters
    folderName = createRunDirectory(nCells,nEdges,nVerts,gamma,lamda,tStar,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,A,B,R)

    # Initialise time and output count
    t = 1E-8
    outputCount = 0
    topologyChange!(A,Ā,Aᵀ,Āᵀ,B,B̄,Bᵀ,B̄ᵀ,C,cellEdgeCount,boundaryVertices)

    while t<tMax

        # 4 step Runge-Kutta integration
        # 1st step of Runge-Kutta
        spatialData!(A,Ā,B,B̄,C,R,nCells,nEdges,cellPositions,cellEdgeCount,cellAreas,cellOrientedAreas,cellPerimeters,cellTensions,cellPressures,edgeLengths,edgeMidpoints,edgeTangents,gamma,preferredPerimeter)
        if t%outputInterval<dt
            visualise(A,Ā,B̄,R,C,F,cellPositions,edgeTangents,edgeMidpoints,nEdges,nVerts,nCells,outputCount,folderName,ϵ,boundaryVertices)
            outputCount+=1
            println("$outputCount/100")
        end
        calculateForce!(F,A,Ā,B,B̄,cellPressures,cellTensions,edgeTangents,edgeLengths,nVerts,nCells,nEdges,ϵ)
        ΔR .= sum(F,dims=2)[:,1,:].*dt/6.0

        # 2nd step of Runge-Kutta
        tempR .= R .+ sum(F,dims=2)[:,1,:].*dt/2.0
        spatialData!(A,Ā,B,B̄,C,tempR,nCells,nEdges,cellPositions,cellEdgeCount,cellAreas,cellOrientedAreas,cellPerimeters,cellTensions,cellPressures,edgeLengths,edgeMidpoints,edgeTangents,gamma,preferredPerimeter)
        calculateForce!(F,A,Ā,B,B̄,cellPressures,cellTensions,edgeTangents,edgeLengths,nVerts,nCells,nEdges,ϵ)
        ΔR .+= sum(F,dims=2)[:,1,:].*dt/3.0

        # 3rd step of Runge-Kutta
        tempR .= R .+ sum(F,dims=2)[:,1,:].*dt/2.0
        spatialData!(A,Ā,B,B̄,C,tempR,nCells,nEdges,cellPositions,cellEdgeCount,cellAreas,cellOrientedAreas,cellPerimeters,cellTensions,cellPressures,edgeLengths,edgeMidpoints,edgeTangents,gamma,preferredPerimeter)
        calculateForce!(F,A,Ā,B,B̄,cellPressures,cellTensions,edgeTangents,edgeLengths,nVerts,nCells,nEdges,ϵ)
        ΔR .+= sum(F,dims=2)[:,1,:].*dt/3.0

        # 4th step of Runge-Kutta
        tempR .= R .+ sum(F,dims=2)[:,1,:].*dt
        spatialData!(A,Ā,B,B̄,C,tempR,nCells,nEdges,cellPositions,cellEdgeCount,cellAreas,cellOrientedAreas,cellPerimeters,cellTensions,cellPressures,edgeLengths,edgeMidpoints,edgeTangents,gamma,preferredPerimeter)
        calculateForce!(F,A,Ā,B,B̄,cellPressures,cellTensions,edgeTangents,edgeLengths,nVerts,nCells,nEdges,ϵ)
        ΔR .+= sum(F,dims=2)[:,1,:].*dt/6.0

        # Result of Runge-Kutta steps
        R .+= ΔR
        t +=dt

    end

    run(`convert -delay 0 -loop 0 output/$folderName/plot"*".png output/$folderName/animated.gif`)

end

export simulate

end
