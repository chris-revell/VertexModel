# Script to minimise the energy of the infinite bilayer model

using Optim
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
using GeometryBasics
using ForwardDiff
using LineSearches

######################################## FUNCTION DEFINITIONS ########################################
function vertexImage(R,κ,λ,s)
    # Function to calculate rotational image of a vector, s=1 for anticlockwise, s=-1 for clockwise
    
    # Set values for images:
    e = 4+κ^2*λ^2
    a = (4-κ^2*λ^2)/e
    b = 4*κ*λ/e
    c = 4λ/e 
    d = 2*κ*λ^2 / e 

    x = R[1] 
    y = R[2]

    R¹ = SVector(x*a -s*y*b + s*c, y*a +s* x*b + d)

    return R¹
end

function repeatedVertexImage(r,κ,λ,n)
    # Apply the periodic image map n times.
    # n > 0: anticlockwise images
    # n < 0: clockwise images
    # n = 0: original vertex

    r_image = r

    if n > 0
        for _ in 1:n
            r_image = vertexImage(r_image, κ, λ, 1)
        end
    elseif n < 0
        for _ in 1:(-n)
            r_image = vertexImage(r_image, κ, λ, -1)
        end
    end

    return r_image
end

function energyCalculation(u,Λ_AA,Λ_BB,Λ_AB)
    # Function to caclulate the energy of the bilayer given input variables x 

    Γ = 0.05

    T = eltype(u)

    # Vertices: [1,2,3,4,5,6,2',3',4'',5''] = [1,2,3,4,5,6,7,8,9,10]
    R = Vector{SVector{2,T}}(undef, 10)
    R[1] = SVector(u[1], u[2])
    R[2] = SVector(u[3], u[4])
    R[3] = SVector(u[5], u[6])
    R[4] = SVector(-u[5], -u[6])   
    R[5] = SVector(u[7], u[8])
    R[6] = SVector(u[9], u[10])
    κ  = u[11]
    λ  = u[12]
    
    ######################### Compute full list of vertices #########################
    R[7] = vertexImage(R[2],κ,λ,1)
    R[8] = vertexImage(R[3],κ,λ,1)
    R[9] = vertexImage(R[4],κ,λ,-1)
    R[10] = vertexImage(R[5],κ,λ,-1)

    ######################### Spatial data for each cell: #########################

    # Order of vertices and edges around each cell (clockewise):
    cellVertexOrders1 = CircularArray([1,7,8,4,3,2])
    cellVertexOrders2 = CircularArray([3,4,5,6,10,9])

    cellEdgeOrders1 = CircularArray([1,2,3,4,5,6])
    cellEdgeOrders2 = CircularArray([4,7,8,9,10,11])

    # Vector of edge tensions for each edge: []
    Λs = [Λ_AA,Λ_AA,Λ_AB,Λ_AB,Λ_AA,Λ_AA,Λ_BB,Λ_BB,Λ_BB,Λ_BB,Λ_AB]

    # Edge lengths
    l_j = Vector{T}(undef, 11)
    l_j[1] = norm(R[7]-R[1])
    l_j[2] = norm(R[8]-R[7])
    l_j[3] = norm(R[4]-R[8])
    l_j[4] = norm(R[4]-R[3])
    l_j[5] = norm(R[2]-R[3])
    l_j[6] = norm(R[1]-R[2])
    l_j[7] = norm(R[5]-R[4])
    l_j[8] = norm(R[6]-R[5])
    l_j[9] = norm(R[10]-R[6])
    l_j[10] = norm(R[9]-R[10])
    l_j[11] = norm(R[3]-R[9])

    # Perimeter of each cell
    L_1 = sum(l_j[cellEdgeOrders1])
    L_2 = sum(l_j[cellEdgeOrders2])
    
    # Area of each cell: 
    A_1 = abs(area(Point{2,T}.(R[cellVertexOrders1])))
    A_2 = abs(area(Point{2,T}.(R[cellVertexOrders2])))

    ######################### Compute energy for each cell #########################
    U = zero(T)
    # Area elasticity term: 
    U += 0.5*((A_1-1)^2 + (A_2-1)^2)
    # Perimeter contractility term: 
    U += 0.5*Γ*(L_1^2+L_2^2)
    # Edge tension term: 
    for j in cellEdgeOrders1
        U += 0.5*Λs[j]*l_j[j]
    end
    for j in cellEdgeOrders2
        U += 0.5*Λs[j]*l_j[j]
    end


    return U
end

function makeBilayerPolygon(R,cellVertexOrder)
    cellPolygons = Vector{Point{2,Float64}}[]
    
    push!(cellPolygons,Point{2,Float64}.(R[cellVertexOrder]))
    
    return cellPolygons
end

function reconstructEquilibrium(x1,y1,x2,y2,x3,y3,x5,y5,x6,y6,κ,λ)
    # Function to reconstruct geometric information at equilibrium 

    # Reconstruct vertices at equilibrium
    R = Vector{SVector{2,Float64}}(undef,10)

    R[1] = SVector(x1,y1)
    R[2] = SVector(x2,y2)
    R[3] = SVector(x3,y3)
    R[4] = SVector(-x3,-y3)
    R[5] = SVector(x5,y5)
    R[6] = SVector(x6,y6)

    R[7]  = vertexImage(R[2],κ,λ,1)
    R[8]  = vertexImage(R[3],κ,λ,1)
    R[9]  = vertexImage(R[4],κ,λ,-1)
    R[10] = vertexImage(R[5],κ,λ,-1)
    
    return R
end

function computeForcesOnInternalVertices(R,Λ_AA,Λ_BB,A,B,κ,λ)

    Γ = 0.05
    Λ_AB = -0.2

    edgeTangents = Vector{SVector{2,Float64}}(undef,19)
    edgeLengths = Vector{Float64}(undef, 19)


    # Find spatial data: 
    edgeTangents[1] = R[7] .- R[1]
    edgeTangents[2] = R[8] .- R[7]
    edgeTangents[3] = R[4] .-  R[8]
    edgeTangents[4] = R[3] .-  R[4]
    edgeTangents[5] = R[2] .-  R[3]
    edgeTangents[6] = R[1] .-  R[2]
    edgeTangents[7] = R[5] .-  R[4]
    edgeTangents[8] = R[6] .-  R[5]
    edgeTangents[9] = R[10] .-  R[6]
    edgeTangents[10] = R[9] .-  R[10]
    edgeTangents[11] = R[3] .-  R[9]
    edgeTangents[12] = vertexImage(R[1],κ,λ,-1) -  vertexImage(R[2],κ,λ,-1)
    edgeTangents[13] = vertexImage(R[2],κ,λ,-1) -  vertexImage(R[3],κ,λ,-1)
    edgeTangents[14] = vertexImage(R[3],κ,λ,-1) -  vertexImage(R[4],κ,λ,-1)
    edgeTangents[15] = vertexImage(R[4],κ,λ,+1) -  vertexImage(R[3],κ,λ,+1)
    edgeTangents[16] = vertexImage(R[5],κ,λ,+1) -  vertexImage(R[4],κ,λ,+1)
    edgeTangents[17] = vertexImage(R[6],κ,λ,+1) -  vertexImage(R[5],κ,λ,+1)
    edgeTangents[18] = vertexImage(R[7],κ,λ,-1) -  vertexImage(R[1],κ,λ,-1)
    edgeTangents[19] = vertexImage(R[10],κ,λ,+1) -  vertexImage(R[6],κ,λ,+1)
    

    # Λs = [Λ_AA, 
    #     Λ_AA,
    #     Λ_AB,
    #     Λ_AB,
    #     Λ_AA,
    #     Λ_AA,
    #     Λ_BB,
    #     Λ_BB,
    #     Λ_BB,
    #     Λ_BB,
    #     Λ_AB
    #     ]

    Λs = zeros(Float64, 19)
    Λs[[1,2,5,6,13,12,18]] .= Λ_AA
    Λs[[7,8,9,10,16,17,19]] .= Λ_BB
    Λs[[3,4,11,14,15]] .= Λ_AB
    
    edgeLengths .= norm.(edgeTangents)

    # Vertices around each cell
    cellVertexOrders1 = CircularArray([1,7,8,4,3,2])
    cellVertexOrders2 = CircularArray([3,4,5,6,10,9])

    cellEdgeOrders1 = CircularArray([1,2,3,4,5,6])
    cellEdgeOrders2 = CircularArray([4,7,8,9,10,11])

    A_1 = abs(area(Point{2,Float64}.(R[cellVertexOrders1])))
    A_2 = abs(area(Point{2,Float64}.(R[cellVertexOrders2])))
    L_1 = sum(edgeLengths[cellEdgeOrders1])
    L_2 = sum(edgeLengths[cellEdgeOrders2])


    cellAreas = [A_1, A_2, A_1, A_2]
    cellPerimeters = [L_1,L_2,L_1,L_2]

    ϵ = SMatrix{2, 2, Float64}([
        0.0 1.0
        -1.0 0.0
    ])

    F = zeros(SVector{2,Float64}, 16, 4)

    Ā = abs.(A)
    B̄ = abs.(B)
    # Writing force on each vertex in terms of incidence matrices instead: 
    for k=1:16
        for j in nzrange(A,k)
            for i in nzrange(B, rowvals(A)[j])
                F[k, rowvals(B)[i]] += 0.5 * (cellAreas[rowvals(B)[i]]-1) * B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ * edgeTangents[rowvals(A)[j]])
                F[k, rowvals(B)[i]] -= Γ*cellPerimeters[rowvals(B)[i]] * B̄[rowvals(B)[i], rowvals(A)[j]] * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
                F[k, rowvals(B)[i]] -= 0.5 * Λs[rowvals(A)[j]] * B̄[rowvals(B)[i], rowvals(A)[j]] *  A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
            end
        end
    end

    
    # Vertex 3: 
    # Cell 1:
    # Edge 4:
    # F[3,1] += 0.5*(A_1 - 1)*(1) .*(ϵ * edgeTangents[4])
    # F[3,1] -= Γ*L_1*(1) .* edgeTangents[4] ./ edgeLengths[4]
    # F[3,1] -= 0.5*Λs[4]*(1) .* edgeTangents[4] ./ edgeLengths[4]
    # # Edge 5: 
    # F[3,1] += 0.5*(A_1 - 1)*(1) .*(ϵ * edgeTangents[5])
    # F[3,1] -= Γ*L_1*(-1) .* edgeTangents[5] ./ edgeLengths[5]
    # F[3,1] -= 0.5*Λs[5]*(-1) .* edgeTangents[5] ./ edgeLengths[5]
    # # Cell 1': 
    # # Edge 5: 
    # F[3,3] += 0.5*(A_1 - 1)*(-1) .*(ϵ * edgeTangents[5])
    # F[3,3] -= Γ*L_1*(-1) .* edgeTangents[5] ./ edgeLengths[5]
    # F[3,3] -= 0.5*Λs[5]*(-1) .* edgeTangents[5] ./ edgeLengths[5]
    # # Edge 11: 
    # F[3,3] += 0.5*(A_1 - 1)*(-1) .*(ϵ * edgeTangents[11])
    # F[3,3] -= Γ*L_1*(1) .* edgeTangents[11] ./ edgeLengths[11]
    # F[3,3] -= 0.5*Λs[11]*(1) .* edgeTangents[11] ./ edgeLengths[11]
    # # Cell 2: 
    # # Edge 11
    # F[3,2] += 0.5*(A_2 - 1)*(1) .*(ϵ * edgeTangents[11])
    # F[3,2] -= Γ*L_2*(1) .* edgeTangents[11] ./ edgeLengths[11]
    # F[3,2] -= 0.5*Λs[11]*(1) .* edgeTangents[11] ./ edgeLengths[11]
    # # Edge 4: 
    # F[3,2] += 0.5*(A_2 - 1)*(-1) .*(ϵ * edgeTangents[4])
    # F[3,2] -= Γ*L_2*(1) .* edgeTangents[4] ./ edgeLengths[4]
    # F[3,2] -= 0.5*Λs[4]*(1) .* edgeTangents[4] ./ edgeLengths[4]

    # # Vertex 4: 
    # # Cell 1:
    # # Edge 3: 
    # F[4,1] += 0.5*(A_1 - 1)*(1) .*(ϵ * edgeTangents[3])
    # F[4,1] -= Γ*L_1*(1) .* edgeTangents[3] ./ edgeLengths[3]
    # F[4,1] -= 0.5*Λs[3]*(1) .* edgeTangents[3] ./ edgeLengths[3]
    # # Edge 4: 
    # F[4,1] += 0.5*(A_1 - 1)*(1) .*(ϵ * edgeTangents[4])
    # F[4,1] -= Γ*L_1*(-1) .* edgeTangents[4] ./ edgeLengths[4]
    # F[4,1] -= 0.5*Λs[4]*(-1) .* edgeTangents[4] ./ edgeLengths[4]
    # # Cell 2: 
    # # Edge 4: 
    # F[4,2] += 0.5*(A_2 - 1)*(-1) .*(ϵ * edgeTangents[4])
    # F[4,2] -= Γ*L_2*(-1) .* edgeTangents[4] ./ edgeLengths[4]
    # F[4,2] -= 0.5*Λs[4]*(-1) .* edgeTangents[4] ./ edgeLengths[4]
    # # Edge 7: 
    # F[4,2] += 0.5*(A_2 - 1)*(1) .*(ϵ * edgeTangents[7])
    # F[4,2] -= Γ*L_2*(-1) .* edgeTangents[7] ./ edgeLengths[7]
    # F[4,2] -= 0.5*Λs[7]*(-1) .* edgeTangents[7] ./ edgeLengths[7]
    # # Cell 2'
    # # Edge 7: 
    # F[4,4] += 0.5*(A_2 - 1)*(-1) .*(ϵ * edgeTangents[7])
    # F[4,4] -= Γ*L_2*(-1) .* edgeTangents[7] ./ edgeLengths[7]
    # F[4,4] -= 0.5*Λs[7]*(-1) .* edgeTangents[7] ./ edgeLengths[7]
    # # Edge 3:
    # F[4,4] += 0.5*(A_2 - 1)*(-1) .*(ϵ * edgeTangents[3])
    # F[4,4] -= Γ*L_2*(1) .* edgeTangents[3] ./ edgeLengths[3]
    # F[4,4] -= 0.5*Λs[3]*(1) .* edgeTangents[3] ./ edgeLengths[3]

    return F
end

function hAroundCellNew!(ii, B, h, ϵ, F, currentNeighbourShell, traversedCells, traversedEdges, cellEdgeOrders, cellVertexOrders)
    # Find an edge belonging to cell ii and shared by any other cell in traversedCells
    if length(traversedCells)>0
        startEdge = (findnz(B[traversedCells, :])[2] ∩ findnz(B[ii, :])[1])[1]
    else 
        startEdge = cellEdgeOrders[ii][1]
    end
    # Find the index of this edge within the ordering of edges around cell ii 
    startInd = findall(x->x==startEdge, cellEdgeOrders[ii])[1]
    # Remember in clockwise ordering, cellEdgeOrders[i][1] precedes cellVertexOrders[1]
    for vertexInd = startInd:(startInd+length(cellVertexOrders[ii])-1)
        h[cellEdgeOrders[ii][vertexInd+1]] = h[cellEdgeOrders[ii][vertexInd]] .+ ϵ*F[cellVertexOrders[ii][vertexInd], ii]
        push!(traversedEdges, cellEdgeOrders[ii][vertexInd]) 
    end
    return nothing 
end

function hNetworkNew(R, A, B, F)
    I = size(B,1)
    J = size(B,2)
    peripheralEdges = findPeripheralEdges(B)
    peripheralCells = findnz(B[:, peripheralEdges.==1])[1]
    cellVertexOrders  = fill(CircularVector(Int64[]), I)
    cellEdgeOrders    = fill(CircularVector(Int64[]), I)
    for i = 1:I
        cellVertexOrders[i], cellEdgeOrders[i] = orderAroundCell(A, B, i)
    end
    ϵ = SMatrix{2, 2, Float64}([
            0.0 1.0
            -1.0 0.0
        ])

    startCell = 1 
    traversedCells = Int64[]
    traversedEdges = Int64[]
    cellNeighbourMatrix = B*transpose(B)
    h = fill(SVector{2, Float64}(zeros(2)), J)

    hAroundCellNew!(startCell, B, h, ϵ, F, Int64[], traversedCells, traversedEdges, cellEdgeOrders, cellVertexOrders)
    push!(traversedCells, startCell)
    while length(traversedCells) < I
        currentNeighbourShell = setdiff(findnz(cellNeighbourMatrix[:, traversedCells])[1], traversedCells)
        for ii in currentNeighbourShell
            hAroundCellNew!(ii, B, h, ϵ, F, currentNeighbourShell, traversedCells, traversedEdges, cellEdgeOrders, cellVertexOrders)
            push!(traversedCells, ii)
        end
    end

    for j=2:J
        h[j] = h[j] - h[1]
    end
    h[1] = @SVector zeros(2)

    return h
end

function curlᵛ(R, A, B, h, cellLinks, E)
    tmp = [-A[j,k]*(cellLinks[j]⋅h[j])/E[k] for j=1:size(A,1), k=[3,4]] # Only at vertices 3 and 4
    return dropdims(sum(tmp, dims=1), dims=1)
end 

function computeCoupleStresses(R,F,κ,λ,A,B)
    # Function to compute couple stresses at interior vertices on the bilayer: 

    h = hNetworkNew(R, A, B, F)
    
    # Initialise couple stress vector for INTERNAL VERTICES ONLY: 
    # coupleStresses = zeros(Float64, 2) # 1: vertex 3, 2: vertex 4

    # CELL CENTRES: Find the centres of cells 1, 2, 3, and 4:
    cellCentres = Vector{SVector{2,Float64}}(undef, 4)
    cellCentres[1] = sum(R[[1,7,8,4,3,2]])/6
    cellCentres[2] = sum(R[[3,4,5,6,10,9]])/6
    cellCentres[3] = vertexImage(cellCentres[1],κ,λ,-1) # periodic image of cell 1 to the left
    cellCentres[4] = vertexImage(cellCentres[2],κ,λ,1)  # periodic image of cell 1 to the right

    # CELL LINKS: only fill in the values relevant to internal vertices 3 and 4.
    cellLinks = zeros(SVector{2,Float64}, 19)
    # Cell links always point from the cell that has B[i,j] = -1 to the cell that has B[i,j]=+1
    cellLinks[5] = cellCentres[1] - cellCentres[3] 
    cellLinks[4] = cellCentres[1] - cellCentres[2]  
    cellLinks[11] = cellCentres[2] - cellCentres[3] 
    cellLinks[3] = cellCentres[1] - cellCentres[4]
    cellLinks[7] = cellCentres[2] - cellCentres[4]   
    
    # println(cellLinks)

    # VERTEX TRIANGLES
    vertexTriangles =  Vector{Point{2,Float64}}[]
    push!(vertexTriangles, Point{2,Float64}.(cellCentres[[1,2,3]])) # vertex 3
    push!(vertexTriangles, Point{2,Float64}.(cellCentres[[1,2,4]])) # vertex 4

    # VERTEX TRIANGLE AREAS: 
    E = zeros(Float64, 16)
    for k=3:4
        E[k] = abs(area(vertexTriangles[k-2]))
    end

    coupleStresses = -curlᵛ(R, A, B, h, cellLinks, E)
    println("coupleStresses = ", coupleStresses)

    return coupleStresses, vertexTriangles

end

function runBilayer(parameterLabelVec,Λ_AAVec,Λ_BBVec,Λ_AB,plotCoupleStresses)

    # flag to plot repeats of the bilayer or not
    plotMultipleImages = true 

    ######################################## INITIALISE VECTORS AND PARAMS ########################################

    # Initialise vectors of each of our minimisation variables 
    x1_vec = []
    y1_vec = []
    x2_vec = []
    y2_vec = []
    x3_vec = []
    y3_vec = []
    x5_vec = []
    y5_vec = []
    x6_vec = []
    y6_vec = []
    κ_vec = []
    λ_vec = []
    γ_AB_vec = []



    # Initialise values for base case 
    λ0 = 1.0
    κ0 = 0.0

    x10 = 1/4
    y10 = 5/(4√3) 
    x20 = -1/4
    y20 = 3/(4√3)
    x30 = -1/4
    y30 = 1/(4√3)
    x50 = 1/4
    y50 = -3/(4√3)
    x60 = -1/4
    y60 = -5/(4√3)

    ########################## Plot initial condition ##########################

    R0 = Vector{SVector{2,Float64}}(undef,10)

    R0[1] = SVector(x10,y10)
    R0[2] = SVector(x20,y20)
    R0[3] = SVector(x30,y30)
    R0[4] = SVector(-x30,-y30)
    R0[5] = SVector(x50,y50)
    R0[6] = SVector(x60,y60)

    R0[7]  = vertexImage(R0[2],κ0,λ0,1)
    R0[8]  = vertexImage(R0[3],κ0,λ0,1)
    R0[9]  = vertexImage(R0[4],κ0,λ0,-1)
    R0[10] = vertexImage(R0[5],κ0,λ0,-1)

    cellVertexOrders1 = CircularArray([1,7,8,4,3,2])
    cellVertexOrders2 = CircularArray([3,4,5,6,10,9])

    polygon1 = makeBilayerPolygon(R0,cellVertexOrders1)
    polygon2 = makeBilayerPolygon(R0,cellVertexOrders2)

    fig_initial = Figure(size=(600,600))
    ax_initial = Axis(fig_initial[1,1], aspect=1)

    hidedecorations!(ax_initial)
    hidespines!(ax_initial)

    poly!(
        ax_initial,
        polygon1,
        color=RGB(255/255,178/255,102/255),
        strokecolor=:black,
        strokewidth=1
    )

    poly!(
        ax_initial,
        polygon2,
        color=RGB(102/255,178/255,255/255),
        strokecolor=:black,
        strokewidth=1
    )

    ax_initial.title = "Initial bilayer geometry"

    save(datadir("InfiniteBilayer","bilayer_initial_geometry.png"), fig_initial)

    # Construct incidence matrices: 

    # Let's try constructing full A and B matrices for the 4-cell system: 
    B = spzeros(Int8, 4, 19)
    A = spzeros(Int8, 19,16)

    B[1,[1,2,3,4,5,6]] .= [1,1,1,1,1,1]
    B[2,[4,7,8,9,10,11]] .= [-1,1,1,1,1,1]
    B[3,[18,5,11,14,13,12]] .= [1,-1,-1,1,1,1]
    B[4,[3,15,16,17,19,7]] .= [-1,1,1,1,1,-1]

    A[1,[1,7]] .= [-1,1]
    A[2,[7,8]] .= [-1,1]
    A[3,[8,4]] .= [-1,1]
    A[4,[4,3]] .= [-1,1]
    A[5,[3,2]] .= [-1,1]
    A[6,[2,1]] .= [-1,1]
    A[7,[4,5]] .= [-1,1]
    A[8,[5,6]] .= [-1,1]
    A[9,[6,10]] .= [-1,1]
    A[10,[10,9]] .= [-1,1]
    A[11,[9,3]] .= [-1,1]
    A[12,[12,11]] .= [-1,1]
    A[13,[13,12]] .= [-1,1]
    A[14,[9,13]] .= [-1,1]
    A[15,[8,14]] .= [-1,1]
    A[16,[14,15]] .= [-1,1]
    A[17,[15,16]] .= [-1,1]
    A[18,[11,2]] .= [-1,1]
    A[19,[16,5]] .= [-1,1]


    for paramInd in 1:length(parameterLabelVec)

        Λ_AA = Λ_AAVec[paramInd]
        Λ_BB = Λ_BBVec[paramInd]

        try

                ######################### INITIALISE BILAYER USING PREVIOUS PARAMS ########################
                u0 = [x10,y10,x20,y20,x30,y30,x50,y50,x60,y60,κ0,λ0]

                # Run minimisation:
                f(u) = energyCalculation(u, Λ_AA, Λ_BB, Λ_AB)

                result = optimize(
                    f,
                    u0,
                    LBFGS();
                    autodiff = :forward,
                )

                u = Optim.minimizer(result)

                x1,y1,x2,y2,x3,y3,x5,y5,x6,y6,κ,λ = u

                push!(x1_vec, x1)
                push!(y1_vec, y1)
                push!(x2_vec, x2)
                push!(y2_vec, y2)
                push!(x3_vec, x3)
                push!(y3_vec, y3)
                push!(x5_vec, x5)
                push!(y5_vec, y5)
                push!(x6_vec, x6)
                push!(y6_vec, y6)
                push!(κ_vec, κ)
                push!(λ_vec, λ)

                # Set initial values for next run: 
                # λ0 = λ
                # κ0 = κ
                # x10 = x1
                # y10 = y1
                # x20 = x2
                # y20 = y2
                # x30 = x3
                # y30 = y3
                # x50 = x5
                # y50 = y5
                # x60 = x6
                # y60 = y6

                push!(γ_AB_vec, Λ_AB - (Λ_AA+Λ_BB)/2)

                ########################## Plot the bilayer at equilibrium ##########################

                # Reconstruct vertices at equilibrium
                R = reconstructEquilibrium(x1, y1, x2, y2, x3, y3, x5, y5, x6, y6, κ, λ)

                F = computeForcesOnInternalVertices(R, Λ_AA, Λ_BB,A,B,κ,λ)

                
                polygon1 = makeBilayerPolygon(R,cellVertexOrders1)
                polygon2 = makeBilayerPolygon(R,cellVertexOrders2)

                if plotMultipleImages
                    fig_bilayer = Figure(size=(900, 700))
                    ax1 = Axis(
                        fig_bilayer[1, 1],
                        aspect = DataAspect()
                    )

                    hidedecorations!(ax1)
                    hidespines!(ax1)

                    # Number of periodic copies on either side of the fundamental unit
                    nRepeats = 10

                    coupleStresses, vertexTriangles = computeCoupleStresses(R,F,κ,λ,A,B)
                    for n in -nRepeats:nRepeats

                        # Apply the periodic transformation to all vertices
                        R_image = [
                            repeatedVertexImage(r, κ, λ, n)
                            for r in R
                        ]

                        polygon1_image = makeBilayerPolygon(
                            R_image,
                            cellVertexOrders1
                        )

                        polygon2_image = makeBilayerPolygon(
                            R_image,
                            cellVertexOrders2
                        )
                        
                        
                        # Highlight the original fundamental unit

                        if n == 0
                            colourA = RGB(255/255, 178/255, 102/255)
                            colourB = RGB(102/255, 178/255, 255/255)
                            lineWidth = 2
                        else
                            colourA = RGBA(255/255, 178/255, 102/255, 0.5)
                            colourB = RGBA(102/255, 178/255, 255/255, 0.5)
                            lineWidth = 1
                        end

                        if plotCoupleStresses
                            vertexTriangleImage1 = [Point{2,Float64}(repeatedVertexImage(p, κ, λ, n)) for p in vertexTriangles[1]]

                            vertexTriangleImage2 = [Point{2,Float64}(repeatedVertexImage(p, κ, λ, n)) for p in vertexTriangles[2]]
                        

                            

                            # Generate a colour map for effective pressures: 
                            cmap3 = cgrad([
                                RGB(0.0, 0.0, 1.0),    # blue
                                RGB(1.0, 1.0, 1.0),   # white, zero
                                RGB(1.0, 0.0, 0.0)   # red
                            ], 256)
                            # Colour bar exclusing exterior vertices 
                            maxabs3 = maximum(abs.(coupleStresses))
                            clims3 = (-maxabs3,maxabs3)
                            # clims3 = (-0.2,0.2)
                            
                            poly!(ax1, vertexTriangleImage1,color = coupleStresses[3], colorrange = clims3, colormap = cmap3,  strokecolor=(:grey, 1.0), strokewidth=0.5)
                            poly!(ax1, vertexTriangleImage2,color = coupleStresses[4], colorrange = clims3, colormap = cmap3,  strokecolor=(:grey, 1.0), strokewidth=0.5)
                            poly!(
                                ax1,
                                polygon1_image,
                                color = :transparent,
                                strokecolor = :black,
                                strokewidth = lineWidth
                            )
                            poly!(
                                ax1,
                                polygon2_image,
                                color = :transparent,
                                strokecolor = :black,
                                strokewidth = lineWidth
                            )

                            # Plot the centre of curvature
                            scatter!(
                                ax1,
                                [0.0],
                                [0.0],
                                color = :red,
                                markersize = 10
                            )

                            Colorbar(fig_bilayer[2,1],colormap = cmap3, colorrange=clims3, label="Couple Stress", width=Relative(0.05),height=Relative(0.6))

                            ax1.title =
                                "Periodic bilayer at the energy minimum for case $(parameterLabelVec[paramInd])"

                            save(
                                datadir(
                                    "InfiniteBilayer", "coupleStressTests",
                                    "bilayer_case_$(parameterLabelVec[paramInd])_periodic_geometry.png"
                                ),
                                fig_bilayer
                            )
                            
                        else
                            poly!(
                                ax1,
                                polygon1_image,
                                color = colourA,
                                strokecolor = :black,
                                strokewidth = lineWidth
                            )

                            poly!(
                                ax1,
                                polygon2_image,
                                color = colourB,
                                strokecolor = :black,
                                strokewidth = lineWidth
                            )

                            # Plot the centre of curvature
                            scatter!(
                                ax1,
                                [0.0],
                                [0.0],
                                color = :red,
                                markersize = 10
                            )
                            ax1.title =
                                "Periodic bilayer at the energy minimum for case $(parameterLabelVec[paramInd])"

                            save(
                                datadir(
                                    "InfiniteBilayer", "multipleRuns",
                                    "bilayer_case_$(parameterLabelVec[paramInd])_periodic_geometry.png"
                                ),
                                fig_bilayer
                            )
                        end

                        
                    end

                else
                    

                    fig_bilayer = Figure(resolution=(600,600))
                    ax1 = Axis(fig_bilayer[1,1],aspect=1)
                    hidedecorations!(ax1)
                    hidespines!(ax1)

                    
                    if plotCoupleStresses 

                        coupleStresses, vertexTriangles = computeCoupleStresses(R,F,κ,λ,A,B)
                        

                        vertexTrianglePoints = Vector{Point{2,Float64}}[]
                        for k =1:2
                            push!(vertexTrianglePoints,Point{2,Float64}.(vertexTriangles[k]))
                        end

                        # Generate a colour map for effective pressures: 
                        cmap3 = cgrad([
                            RGB(0.0, 0.0, 1.0),    # blue
                            RGB(1.0, 1.0, 1.0),   # white, zero
                            RGB(1.0, 0.0, 0.0)   # red
                        ], 256)
                        # Colour bar exclusing exterior vertices 
                        # maxabs3 = maximum(abs.(coupleStresses))
                        # clims3 = (-maxabs3,maxabs3)
                        clims3 = (-0.1,0.1)

                        for k =1:2
                            poly!(ax1, vertexTrianglePoints[k],color = coupleStresses[k], colorrange = clims3, colormap = cmap3,  strokecolor=(:grey, 1.0), strokewidth=0.5)
                        end

                        poly!(ax1, polygon1, color=:transparent, strokecolor=(:black, 1.0), strokewidth=1)
                        poly!(ax1, polygon2, color=:transparent, strokecolor=(:black, 1.0), strokewidth=1)

                    else
                        poly!(ax1, polygon1, color=RGB(255/255,178/255,102/255), strokecolor=(:black, 1.0), strokewidth=1)
                        poly!(ax1, polygon2, color=RGB(102/255,178/255,255/255), strokecolor=(:black, 1.0), strokewidth=1)
                    end

                    
                    scatter!(
                        ax1,
                        [r[1] for r in R],
                        [r[2] for r in R],
                        color=:black
                    )
                    scatter!(
                        ax1,
                        [0],
                        [0],
                        color=:red
                    )
                    text!(
                        ax1,
                        [0],
                        [0],
                        text="O",
                        position=:lt,
                        align=(:left,:top),
                        color=:red
                    )
                    text!(
                        ax1,
                        [r[1] for r in R],
                        [r[2] for r in R],
                        text=string.(1:10)
                    )

                    ax1.title = "Bilayer geometry at energy minimum for case $(parameterLabelVec[paramInd])"
                    save(datadir("InfiniteBilayer","multipleRuns","bilayer_case_$(parameterLabelVec[paramInd])_geometry_coupleStress.png"), fig_bilayer)

                end
        
        catch err

                    γ_AB = Λ_AB - (Λ_AA + Λ_BB) / 2

                    # Preserve vector lengths and parameter indexing
                    push!(x1_vec, NaN)
                    push!(y1_vec, NaN)
                    push!(x2_vec, NaN)
                    push!(y2_vec, NaN)
                    push!(x3_vec, NaN)
                    push!(y3_vec, NaN)
                    push!(x5_vec, NaN)
                    push!(y5_vec, NaN)
                    push!(x6_vec, NaN)
                    push!(y6_vec, NaN)
                    push!(κ_vec, NaN)
                    push!(λ_vec, NaN)
                    push!(γ_AB_vec, γ_AB)
            
                    println(
                        stderr,
                        "Failed $paramInd/$(length(parameterLabelVec)): ",
                        "Λ_AA = $Λ_AA, Λ_BB = $Λ_BB"
                    )
                    showerror(stderr, err)
                    println(stderr)    



        
        end
        
        
        
    end

    return x1_vec, y1_vec, x2_vec, y2_vec, x3_vec, y3_vec, x5_vec, y5_vec, x6_vec, y6_vec, κ_vec, λ_vec,γ_AB_vec
end



# Λ_AAVec = [-0.3389236545682103]
# Λ_BBVec = [-0.10480831475885634]
Λ_AAVec = [-0.2,-0.1,-0.1,-0.15,-0.3,-0.3,-0.3,-0.25,
# -0.339,
# -0.105,
-0.2,-0.2,-0.25,-0.339,
-0.3,
-0.3
]
Λ_BBVec = [-0.2,-0.25,-0.3,-0.3,-0.3,-0.15,-0.1,-0.1,
# -0.105,
# -0.339,
-0.339,
-0.3,
-0.3,
-0.2,-0.2,-0.25
]
# Λ_AAVec = [-0.01]
# Λ_BBVec = [-0.38]
Λ_AB = -0.2

# parameterLabelVec = ["(O)","(P)","(Q)","(R)","(S)","test","test2"]
parameterLabelVec = []
for (ind,AAval) in enumerate(Λ_AAVec)
    BBval = Λ_BBVec[ind]
    push!(parameterLabelVec, ["(Λ_{AA},Λ_{BB}) = ($AAval,$BBval)"])
end

plotCoupleStresses=true

x1_vec, y1_vec, x2_vec, y2_vec, x3_vec, y3_vec, x5_vec, y5_vec, x6_vec, y6_vec, κ_vec, λ_vec,γ_AB_vec = runBilayer(parameterLabelVec,Λ_AAVec,Λ_BBVec,Λ_AB,plotCoupleStresses)


