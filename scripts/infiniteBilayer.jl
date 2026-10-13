# Code to set up the infinite bilayer problem. 

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


include(srcdir("Callbacks.jl")); using .Callbacks 

# PICK PARAMETERS 


# INITIALISE MONOLAYER 
    #   - Pick initial distances of points, then define where we want the periodic boundary to be
    #   - Define the initial state vector to be the points that we can see inside the monolayer 
    #   - Need to use kappa and lambda when we are working out edge vectors and stuff, the same way we did on the periodc box case

    #   - Calculations of spatial data will involve our parameters lambda and kappa. 

    # R = vector of vertex positions, 1:6
    # u0 = initial guess of minimisation variables [x₁,y₁,x₂,y₂,x₃,y₃,x₅,y₅,x₆,y₆,κ,λ]


Δy =sqrt(3)/2
R = SVector{2,Float64}[]

# Top cell: 
push!(R,SVector(0.0, Δy))   #1
push!(R,SVector(-0.5, Δy/2))  #2
push!(R, SVector(-0.5, -Δy/2))    # 3
# Interface
push!(R, SVector( 0.5, -Δy/2))    # 4
# Bottom cell
push!(R, SVector( 0.5, -3Δy/2))   # 5
push!(R, SVector( 0.0, -2Δy))     # 6

params.nVerts = 6
params.nEdges = 8
params.nCells = 2

κ0 = 0
λ0 = 1.0 # Cell width for how we have initialised the monolayer 

# Initialise edges in incidence matrices: 
A = spzeros(Int64, nEdges, nVerts)
B = spzeros(Int64, nCells, nEdges)

A[1, 1] = -1; A[1, 2] = 1
A[2, 2] = -1; A[2, 3] = 1
A[3, 3] = -1; A[3, 4] = 1
A[4, 4] = -1; A[4, 3] = 1
A[5, 4] = -1; A[5, 5] = 1
A[6, 5] = -1; A[6, 6] = 1
A[7, 6] = -1; A[7, 5] = 1
A[8, 1] = -1; A[8, 2] = 1

B[1, 1] = -1; B[1, 2] = -1; B[1, 3] = -1; B[1, 4] = -1; B[1, 8] = 1
B[2, 3] = 1; B[2, 4] = 1; B[2, 5] = 1; B[2, 6] = 1; B[2, 7] = 1

# Start by assigning all cells as type A:
cellsTypeB = [2]
cellsTypeA = [1]
# Label A cells as 0, B cells as 1.
cellLabels = zeros(Int64, nCells)
cellLabels[cellsTypeB] .= 1
cellL₀s = L₀.*ones(nCells)

# Fill preallocated matrices into struct for convenience
matrices = MatricesContainer(
    A                 = A,
    B                 = B,
    Aᵀ                = spzeros(Int64, nVerts, nEdges),
    Ā                 = spzeros(Int64, nEdges, nVerts),
    Āᵀ                = spzeros(Int64, nVerts, nEdges),
    Bᵀ                = spzeros(Int64, nEdges, nCells),
    B̄                 = spzeros(Int64, nCells, nEdges),
    B̄ᵀ                = spzeros(Int64, nEdges, nCells),
    C                 = spzeros(Int64, nCells, nVerts),
    cellEdgeCount     = zeros(Int64, nCells),
    cellVertexOrders  = fill(CircularVector(Int64[]), nCells),
    cellEdgeOrders    = fill(CircularVector(Int64[]), nCells),
    boundaryVertices  = zeros(Int64, nVerts),
    boundaryEdges     = zeros(Int64, nEdges),
    boundaryCells     = zeros(Int64, nCells),
    cellPositions     = fill(SVector{2,Float64}(zeros(2)), nCells),
    cellPerimeters    = zeros(nCells),
    cellOrientedAreas = fill(SMatrix{2,2,Float64}(zeros(2,2)), nCells),
    cellAreas         = zeros(nCells),
    cellA₀s           = fill(A₀, nCells),
    cellL₀s           = cellL₀s,
    cellTensions      = zeros(nCells),
    cellPressures     = zeros(nCells),
    cellTimeToDivide  = cellTimeToDivide,
    μ                 = ones(nCells),
    Γ                 = γ.*ones(nCells),
    edgeLengths       = zeros(nEdges),
    edgeTangents      = fill(SVector{2,Float64}(zeros(2)), nEdges),
    edgeMidpoints     = fill(SVector{2,Float64}(zeros(2)), nEdges),
    edgeMidpointLinks = spzeros(SVector{2,Float64}, nCells, nVerts),
    timeSinceT1       = zeros(nEdges),
    vertexAreas       = ones(nVerts),
    F                 = spzeros(SVector{2,Float64}, nVerts, nCells),
    externalF         = fill(SVector{2,Float64}(zeros(2)), nVerts),
    totalF            = fill(SVector{2,Float64}(zeros(2)), nVerts),
    ϵ                 = SMatrix{2, 2, Float64}([
                            0.0 1.0
                            -1.0 0.0
                        ]), # clockwise orientation
    cellShapeTensor   = fill(SMatrix{2,2,Float64}(zeros(2,2)), nCells),
    cellLabels        = cellLabels,
    Λs                = zeros(nEdges),
    firstT1onEdge     = zeros(Int64, nEdges),
    P_effs            = zeros(nCells),
    T_effs            = zeros(nCells),
    ξs                = zeros(nCells),
    e₁                = fill(SVector{2,Float64}(zeros(2)), nCells),
    edgeLabels        = zeros(Int64, nEdges),
    trackedVertDistance = [],
    trackedTimePoints = [],
    jsAfterAblation   = [],
)

# Pack parameters into a struct for convenience
params = ParametersContainer(
    initialSystem     = initialSystem,
    boundaryType      = boundaryType,
    cellLayout        = cellLayout,
    nCells            = nCells,
    nEdges            = nEdges,
    nVerts            = nVerts,
    γ                 = γ,
    L₀                = L₀,
    A₀                = A₀,
    pressureExternal  = pressureExternal,
    outputTotal       = outputTotal,
    outputInterval    = outputInterval,
    viscousTimeScale  = viscousTimeScale,
    realTimetMax      = realTimetMax,
    tMax              = tMax,
    realCycleTime     = realCycleTime,
    nCycles           = nCycles,
    nonDimCycleTime   = nonDimCycleTime,
    t1Threshold       = 0,
    peripheralTension = peripheralTension,
    β                 = β,
    seed              = seed,
    rng               = rng,
    distLogNormal     = LogNormal(0.0, 0.2),
    energyModel       = energyModel,
    vertexWeighting   = vertexWeighting,
    noiseWeighting    = noiseWeighting,
    cellsTypeA        = cellsTypeA,
    cellsTypeB        = cellsTypeB,
    L_x               = L_x,
    L_y               = L_y,
    Λ_AA              = Λ_AA,
    Λ_AB              = Λ_AB,
    Λ_BB              = Λ_BB,
    # Λ_AE              = 0.5*Λ_AA,
    # Λ_BE              = 0.5*Λ_BB,
    Λ_AE              = Λ_AE,
    Λ_BE              = Λ_BE,
    Area_A_ratio      = Area_A_ratio,
    t1timeGap         = t1timeGap,
    l_AA              = l_AA,
    l_BB              = l_BB,
    l_AB              = l_AB,
    l_AE              = l_AE,
    l_BE              = l_BE,
    k_tracked         = zeros(2),
    jAblated         = 0,
)

# Initial evaluation of matrices based on system topology
topologyChange!(R,params,matrices)
spatialDataBilayer!(R, params, matrices)

u0 = float64[] # Initial guess vector 
for (rInd,r) in enumerate(R) 
    if rInd != 4 # We skip this since it is defined by point 3 and the origin
        push!(u0, r[1])
        push!(u0, r[2])
    end
end
# Add κ and λ
push!(u0, κ0)
push!(u0, λ0) 


abstol = 1e-9
reltol = 1e-6
sstol = 100.0*abstol # Steady state tolerance for callback
solver = Tsit5()
# Set ODE parameters: 
prob = ODEProblem(bilayerModel!,
        u0,
        (0.0, (termSteadyState ? Inf : params.tMax)),
        (termSteadyState ? (params, matrices, sstol) : (params, matrices)),
    )
alltStops = collect(0.0:params.outputInterval:params.tMax)# Time points beyond which we plot the monolayer
integrator = init(prob,
    solver,
    tstops=alltStops,
    abstol=abstol,
    reltol=reltol,
    save_on=false,
    save_start=false,
    save_end=true,
    callback = DiscreteCallback(termSteadyState ? conditionSteadyState : conditiontMax, terminate!),
        )  

while (integrator.sol.retcode == ReturnCode.Default || integrator.sol.retcode == ReturnCode.Success) && integrator.sol.retcode !=:Terminate


    optVec = reinterpret(SVector{Float64},integrator.u)

    if any(!isfinite, integrator.u)
        @show integrator.t
        @show integrator.u
        error("NaN or Inf detected in integrator.u")
    end

    # Step integrator forwards in time to update minimisation variables (vertex positions, κ, λ)
    step!(integrator)

    
 

    # Update spatial data (edge lengths, cell areas, etc.) following iteration of the integrator
    spatialDataBilayer!(R, params, matrices)
        # This is where we will update the value of lambda and kappa based on how the vertices have changed
        # Is that right? Let's try.
end
