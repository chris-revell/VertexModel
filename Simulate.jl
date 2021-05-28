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
using InitialHexagons
using Visualise
using T1Transitions
using Division

@inline @views function simulate(initialSystem)

    # Parameters
    realTimetMax       = 10000.0 # Real time maximum system run time /seconds
    gamma              = 0.2    # Parameters in energy relaxation
    lamda              = -0.3   # Parameters in energy relaxation
    tStar              = 20.0   # Relaxation rate, approx from Sarah's data.
    dt                 = 0.01   # Non dimensionalised time step
    preferredArea      = 1.0    # Cell preferred area (1.0 by default)
    pressureExternal   = 0.1    # External pressure applied isotropically to system boundary
    cellCycleTime      = 4000.0 # Cell cycle time in seconds
    outputTotal        = 20                  # Number of data outputs
    t1Threshold        = 0.01   # Edge length at which a T1 transition is triggered

    # Derived parameters
    tMax               = realTimetMax/tStar  # Non dimensionalised maximum system run time
    outputInterval     = tMax/outputTotal    # Time interval for storing system data (non dimensionalised)
    preferredPerimeter = -lamda/(2*gamma)    # Cell preferred perimeter
    nonDimCycleTime    = cellCycleTime/tStar # Non dimensionalised cell cycle time
    ϵ                  = [0.0 1.0
                         -1.0 0.0]           # Antisymmetric rotation matrix

    # Initialise system matrices from function or file
    # A = Incidence matrix. Rows => edges; columns => vertices.
    # B = Incidence matrix. Rows => cells; columns => edges. Values +/-1 for orientation
    # (Derived) C = adjacency matrix. Rows => cells; Columns => vertices. = 0.5*B̄*Ā
    # R = Coordinates of vertices
    if initialSystem=="single"
        A,B,R = initialHexagons(1)
    elseif initialSystem=="triple"
        A,B,R = initialHexagons(3)
    else
        # Import system matrices from file
        A = sparse(readdlm("input/$(initialSystem)_A.txt",',',Int64,'\n'))
        B = sparse(readdlm("input/$(initialSystem)_B.txt",',',Int64,'\n'))
        R = readdlm("input/$(initialSystem)_R.txt",',',Float64,'\n')
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
    cellAges          = rand(nCells).*nonDimCycleTime# 1D vector of initial cell ages.
    edgeLengths       = zeros(nEdges,1)              # 1D matrix of scalar edge lengths
    edgeTangents      = zeros(nEdges,2)              # 2D matrix of tangent vectors for each edge (magnitude = edge length)
    edgeMidpoints     = zeros(nEdges,2)              # 2D matrix of position coordinates for each edge midpoint
    vertexEdges       = zeros(Int64,nVerts,3)        # 2D matrix containing the labels of all 3 edges around each vertex
    vertexCells       = zeros(Int64,nVerts,3)        # 2D matrix containing the labels of all 2-3 cells around each vertex
    F                 = zeros(nVerts,2)              # 3D array containing force vectors on vertex k from cell i, Fᵢₖ
    Fexternal         = zeros(nVerts,2)              # 3D array containing force vectors on vertex k from cell i, Fᵢₖ

    # Create output directory in which to store results and parameters
    folderName = createRunDirectory(nCells,nEdges,nVerts,gamma,lamda,tStar,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,A,B,R)

    # Initialise time and output count
    t = 1E-8   # Initial time is very small but slightly above 0 to avoid issues with remainders in output interval calculation
    outputCount = 0
    transitionOccurred = 0

    topologyChange!(A,Ā,Aᵀ,Āᵀ,B,B̄,Bᵀ,B̄ᵀ,C,R,cellEdgeCount,cellPositions,boundaryVertices,vertexEdges,edgeTangents,nVerts,nCells)

    while t<tMax

        if maximum(cellAges) > nonDimCycleTime
            A, B, nCells, nEdges = division!(A,B,C,R,cellAges,cellEdgeCount,cellPositions,edgeMidpoints,nonDimCycleTime,nCells,nEdges,nVerts)
            append!(cellAges, zeros(nCells-length(cellEdgeCount)))
            Aᵀ                = spzeros(Int64,nVerts,nEdges)
            Ā                 = spzeros(Int64,nEdges,nVerts)
            Āᵀ                = spzeros(Int64,nVerts,nEdges)
            Bᵀ                = spzeros(Int64,nEdges,nCells)
            B̄                 = spzeros(Int64,nCells,nEdges)
            B̄ᵀ                = spzeros(Int64,nEdges,nCells)
            C                 = spzeros(Int64,nCells,nVerts)
            tempR             = zeros(nVerts,2)
            ΔR                = zeros(nVerts,2)
            cellEdgeCount     = zeros(Int64,nCells,1)
            boundaryVertices  = zeros(Int64,nVerts,1)
            cellPositions     = zeros(nCells,2)
            cellPerimeters    = zeros(nCells,1)
            cellOrientedAreas = zeros(nCells,2,2)
            cellAreas         = zeros(nCells,1)
            cellTensions      = zeros(nCells,1)
            cellPressures     = zeros(nCells,1)
            edgeLengths       = zeros(nEdges,1)
            edgeTangents      = zeros(nEdges,2)
            edgeMidpoints     = zeros(nEdges,2)
            vertexEdges       = zeros(Int64,nVerts,3)
            vertexCells       = zeros(Int64,nVerts,3)
            F                 = zeros(nVerts,2)
            Fexternal         = zeros(nVerts,2)
            topologyChange!(A,Ā,Aᵀ,Āᵀ,B,B̄,Bᵀ,B̄ᵀ,C,R,cellEdgeCount,cellPositions,boundaryVertices,vertexEdges,edgeTangents,nVerts,nCells)
        end

        # 4 step Runge-Kutta integration
        # 1st step of Runge-Kutta
        spatialData!(A,Ā,B,B̄,C,R,nCells,nEdges,cellPositions,cellEdgeCount,cellAreas,cellOrientedAreas,cellPerimeters,cellTensions,cellPressures,edgeLengths,edgeMidpoints,edgeTangents,gamma,preferredPerimeter)
        transitionOccurred = t1Transitions!(A,Ā,B,B̄,C,R,nEdges,edgeLengths,edgeTangents,t1Threshold,ϵ)
        if transitionOccurred==1
            topologyChange!(A,Ā,Aᵀ,Āᵀ,B,B̄,Bᵀ,B̄ᵀ,C,R,cellEdgeCount,cellPositions,boundaryVertices,vertexEdges,edgeTangents,nVerts,nCells)
            spatialData!(A,Ā,B,B̄,C,R,nCells,nEdges,cellPositions,cellEdgeCount,cellAreas,cellOrientedAreas,cellPerimeters,cellTensions,cellPressures,edgeLengths,edgeMidpoints,edgeTangents,gamma,preferredPerimeter)
            transitionOccurred = 0
        end
        if t%outputInterval<dt
            visualise(A,Ā,B̄,R,C,F,cellPositions,edgeTangents,edgeMidpoints,nEdges,nVerts,nCells,outputCount,folderName,ϵ,boundaryVertices,vertexEdges)
            outputCount+=1
            println("$outputCount/$outputTotal")
        end
        calculateForce!(F,Fexternal,A,Ā,B,B̄,cellPressures,cellTensions,edgeTangents,edgeLengths,nVerts,nCells,nEdges,ϵ,pressureExternal,boundaryVertices)
        ΔR .= (F .+ Fexternal).*dt/6.0

        # 2nd step of Runge-Kutta
        tempR .= R .+ (F .+ Fexternal).*dt/2.0
        spatialData!(A,Ā,B,B̄,C,tempR,nCells,nEdges,cellPositions,cellEdgeCount,cellAreas,cellOrientedAreas,cellPerimeters,cellTensions,cellPressures,edgeLengths,edgeMidpoints,edgeTangents,gamma,preferredPerimeter)
        calculateForce!(F,Fexternal,A,Ā,B,B̄,cellPressures,cellTensions,edgeTangents,edgeLengths,nVerts,nCells,nEdges,ϵ,pressureExternal,boundaryVertices)
        ΔR .+= (F .+ Fexternal).*dt/3.0

        # 3rd step of Runge-Kutta
        tempR .= R .+ (F .+ Fexternal).*dt/2.0
        spatialData!(A,Ā,B,B̄,C,tempR,nCells,nEdges,cellPositions,cellEdgeCount,cellAreas,cellOrientedAreas,cellPerimeters,cellTensions,cellPressures,edgeLengths,edgeMidpoints,edgeTangents,gamma,preferredPerimeter)
        calculateForce!(F,Fexternal,A,Ā,B,B̄,cellPressures,cellTensions,edgeTangents,edgeLengths,nVerts,nCells,nEdges,ϵ,pressureExternal,boundaryVertices)
        ΔR .+= (F .+ Fexternal).*dt/3.0

        # 4th step of Runge-Kutta
        tempR .= R .+ (F .+ Fexternal).*dt
        spatialData!(A,Ā,B,B̄,C,tempR,nCells,nEdges,cellPositions,cellEdgeCount,cellAreas,cellOrientedAreas,cellPerimeters,cellTensions,cellPressures,edgeLengths,edgeMidpoints,edgeTangents,gamma,preferredPerimeter)
        calculateForce!(F,Fexternal,A,Ā,B,B̄,cellPressures,cellTensions,edgeTangents,edgeLengths,nVerts,nCells,nEdges,ϵ,pressureExternal,boundaryVertices)
        ΔR .+= (F .+ Fexternal).*dt/6.0

        # Result of Runge-Kutta steps
        R .+= ΔR
        t +=dt
        cellAges .+= dt

    end

    run(`convert -delay 0 -loop 0 output/$folderName/plot"*".png output/$folderName/animated.gif`)

end

export simulate

end
