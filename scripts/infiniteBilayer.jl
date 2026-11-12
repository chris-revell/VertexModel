# Script to minimise the energy of the infinite bilayer model and compute couple stresses at the composite interface

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

include(srcdir("PlotSetup.jl")); using .PlotSetup
include(srcdir("Visualise.jl")); using .Visualise

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
    cellVertexOrders2 = CircularArray([4,5,6,10,9,3])

    cellEdgeOrders1 = CircularArray([6,1,2,3,4,5])
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

function constructIncidenceMatrices()
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

    return A,B
end

function computeForces(R_ext,Λ_AA,Λ_BB,Λ_AB,Γ,A,B)
    edgeTangents = findEdgeTangents(R_ext,A)
    edgeLengths = findEdgeLengths(R_ext,A)

    Λs = zeros(Float64, 19)
    Λs[[1,2,5,6,13,12,18]] .= Λ_AA
    Λs[[7,8,9,10,16,17,19]] .= Λ_BB
    Λs[[3,4,11,14,15]] .= Λ_AB

    cellAreas = findCellAreas(R_ext,A,B)
    cellPerimeters = findCellPerimeterLengths(R_ext,A,B)

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
    return F
end

function computeCoupleStresses(R,κ,λ,A,B,Λ_AA,Λ_BB,Λ_AB,Γ)

    # Construct the full list of vertices for the 4-cell system: 
    R_ext = Vector{SVector{2,Float64}}(undef,16)
    R_ext[1:10] = R 
    R_ext[11] = vertexImage(R[1],κ,λ,-1)
    R_ext[12] = vertexImage(R[2],κ,λ,-1)
    R_ext[13] = vertexImage(R[3],κ,λ,-1)
    R_ext[14] = vertexImage(R[4],κ,λ,+1)
    R_ext[15] = vertexImage(R[5],κ,λ,+1)
    R_ext[16] = vertexImage(R[6],κ,λ,+1)

    # println(R_ext)

    F = computeForces(R_ext,Λ_AA,Λ_BB,Λ_AB,Γ,A,B)



    h = hNetwork(R_ext,A,B,F)

    coupleStresses = -curlᵛ(R_ext,A,B,h)

    vertexTriangles = findCellLinkVertexTriangles(R_ext,A,B)

    return coupleStresses, vertexTriangles, R_ext
end

function checkHNetwork(R, A, B, F, Λ_AA,Λ_BB)
    # Function to work with DiscreteCalculus.jl couple stress functions, from a simulation in equilbrium. 
    

    h = hNetwork(R,A,B,F)
    coupleStresses = -curlᵛ(R,A,B,h)
    vertexTriangles = findCellLinkVertexTriangles(R,A,B)

    cellEdgeOrders = Dict{Int64,CircularArray{Int64}}()
    cellVertexOrders = Dict{Int64,CircularArray{Int64}}()

    cellEdgeOrders[1] = CircularArray([6,1,2,3,4,5])
    cellVertexOrders[1] = CircularArray([1,7,8,4,3,2])

    cellEdgeOrders[2] = CircularArray([4,7,8,9,10,11])
    cellVertexOrders[2] = CircularArray([4,5,6,10,9,3])

    cellEdgeOrders[3] = CircularArray([5,11,14,13,12,18])
    cellVertexOrders[3] = CircularArray([3,9,13,12,11,2])

    cellEdgeOrders[4] = CircularArray([3,15,16,17,19,7])
    cellVertexOrders[4] = CircularArray([8,14,15,16,5,4])
    

    ϵ = SMatrix{2, 2, Float64}([
        0.0 1.0
        -1.0 0.0
    ]) 
    
    # Case of plotting a cluster of cells to check the rotated forces align with the vector force potential; 
    forceComparisonFig, forceAx, polyAx, forceComparisonCbar = forceComparisonPlotSetup()

    scatter!(
        forceAx,
        Point2f.(h),
        color = :blue
    )

    for j in eachindex(h)
        scatter!(
            forceAx,
            [Point2f(h[j])],
            color = :blue,
            markersize = 12
        )
    
        text!(
            forceAx,
            string(j),
            position = Point2f(h[j]),
            offset = (8, 8),
            color = :black,
            fontsize = 14
        )
    end

    clusterCells = [1,2,3,4]
    # clusterCells=[3]
    startCell = clusterCells[1]
    traversedCells = Int64[]
    traversedEdges = Int64[]
    cellNeighbourMatrix = B*transpose(B)
    rotatedVerticesList = fill(SVector{2, Float64}(zeros(2)), size(B,2)) # This is what we will compare with h
    
    # cellColours = Dict(zip(clusterCells, distinguishable_colors(length(clusterCells))))

    startEdge = cellEdgeOrders[startCell][1]
    # Find the index of this edge within the ordering of edges around cell ii 
    startInd = findall(x->x==startEdge, cellEdgeOrders[startCell])[1]
    for vertexInd = startInd:(startInd+length(cellVertexOrders[startCell])-1)
        if length(traversedEdges)==0
            rotatedVerticesList[cellEdgeOrders[startCell][startInd]] = h[cellEdgeOrders[startCell][startInd]]
        end
        rotatedVerticesList[cellEdgeOrders[startCell][vertexInd+1]] = rotatedVerticesList[cellEdgeOrders[startCell][vertexInd]] .+ ϵ*F[cellVertexOrders[startCell][vertexInd], startCell]

        arrows!(forceAx,[Point2f(rotatedVerticesList[cellEdgeOrders[startCell][vertexInd]])],[Vec2f(ϵ*F[cellVertexOrders[startCell][vertexInd], startCell])],color=:red, linewidth=1, arrowsize=0.03)
        
        push!(traversedEdges, cellEdgeOrders[startCell][vertexInd]) 
    end
    push!(traversedCells, startCell)

    while length(traversedCells) < length(clusterCells)
        currentNeighbourShell = setdiff(findnz(cellNeighbourMatrix[:, traversedCells])[1], traversedCells)
        for ii in currentNeighbourShell 
            if ii in clusterCells

                startEdge = cellEdgeOrders[ii][1]
                # Find the index of this edge within the ordering of edges around cell ii 
                startInd = findall(x->x==startEdge, cellEdgeOrders[ii])[1]

                rotatedVerticesList[cellEdgeOrders[ii][startInd]] = h[cellEdgeOrders[ii][startInd]]
                # Remember in clockwise ordering, cellEdgeOrders[i][1] precedes cellVertexOrders[1]
                for vertexInd = startInd:(startInd+length(cellVertexOrders[ii])-1)
                    rotatedVerticesList[cellEdgeOrders[ii][vertexInd+1]] = rotatedVerticesList[cellEdgeOrders[ii][vertexInd]] .+ ϵ*F[cellVertexOrders[ii][vertexInd], ii]
                    push!(traversedEdges, cellEdgeOrders[ii][vertexInd]) 
                    arrows!(forceAx,[Point2f(rotatedVerticesList[cellEdgeOrders[ii][vertexInd]])],[Vec2f(ϵ*F[cellVertexOrders[ii][vertexInd], ii])],color=:red, linewidth=1, arrowsize=0.03)
                end
                
                push!(traversedCells, ii)
            end
        end
    end

    # forceComparisonCbar = visualiseForceComparison(R, h, clusterCells, forceComparisonFig, polyAx,forceAx, forceComparisonCbar, params, matrices, coupleStresses, vertexTriangles)
    save(datadir("InfiniteBilayer","coupleStressTests", "forceComparisonPlot($Λ_AA,$Λ_BB).png"), forceComparisonFig)


    return nothing
end

##################################### Main funciton ###############################################

function runBilayer(Λ_vec, Λ_AB, Γ, nRepeats)

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

    # Construct incidence matrices
    A,B = constructIncidenceMatrices()
    for paramInd in eachindex(Λ_vec)
        Λ_AA = Λ_vec[paramInd][1]
        Λ_BB = Λ_vec[paramInd][2]

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
            
            
    
            ########################## Plot the bilayer at equilibrium ##########################
    
            # Reconstruct vertices at equilibrium
            R = reconstructEquilibrium(x1, y1, x2, y2, x3, y3, x5, y5, x6, y6, κ, λ)

            # Order of vertices and edges around each cell (clockewise):
            cellVertexOrders1 = CircularArray([1,7,8,4,3,2])
            cellVertexOrders2 = CircularArray([4,5,6,10,9,3])
    
            polygon1 = makeBilayerPolygon(R,cellVertexOrders1)
            polygon2 = makeBilayerPolygon(R,cellVertexOrders2)
    
            fig_bilayer = Figure(size=(900, 700))
            ax1 = Axis(
                fig_bilayer[1, 1],
                aspect = DataAspect()
            )
            ax2 = Axis(
                fig_bilayer[2, 1],
                aspect = DataAspect()
            )
    
            hidedecorations!(ax1)
            hidespines!(ax1)
            hidedecorations!(ax2)
            hidespines!(ax2)
    
            # Compute couple stresses: 
            coupleStresses, vertexTriangles, R_ext = computeCoupleStresses(R,κ,λ,A,B,Λ_AA,Λ_BB,Λ_AB,Γ)

            # println(abs.(area.(vertexTriangles)))

            # Generate a colour map for effective pressures: 
            cmap3 = cgrad([
                RGB(0.0, 0.0, 1.0),    # blue
                RGB(1.0, 1.0, 1.0),   # white, zero
                RGB(1.0, 0.0, 0.0)   # red
            ], 256)
            # Colour bar excluding exterior vertices 
            maxabs3 = maximum(abs.(coupleStresses[[3,4]]))
            clims3 = (-maxabs3,maxabs3)
            # clims3 = (-0.02,0.02)
    
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
    
                ########################### PLOT POLYGONS BY CELL COLOUR ##########################
    
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
    
                ########################## COUPLE STRESS PLOT: ##########################
    
                vertexTriangleImage3 = [Point{2,Float64}(repeatedVertexImage(p, κ, λ, n)) for p in vertexTriangles[3]]
    
                vertexTriangleImage4 = [Point{2,Float64}(repeatedVertexImage(p, κ, λ, n)) for p in vertexTriangles[4]]
    
                poly!(ax2, vertexTriangleImage3,color = coupleStresses[3], colorrange = clims3, colormap = cmap3,  strokecolor=(:grey, 1.0), strokewidth=0.5)
                poly!(ax2, vertexTriangleImage4,color = coupleStresses[4], colorrange = clims3, colormap = cmap3,  strokecolor=(:grey, 1.0), strokewidth=0.5)
                
                poly!(
                    ax2,
                    polygon1_image,
                    color = :transparent,
                    strokecolor = :black,
                    strokewidth = lineWidth
                )
                poly!(
                    ax2,
                    polygon2_image,
                    color = :transparent,
                    strokecolor = :black,
                    strokewidth = lineWidth
                )
    
                
            end
    
            # Plot the centre of curvature
            scatter!(
                ax1,
                [0.0],
                [0.0],
                color = :red,
                markersize = 10
            )
            ax1.title =
                "Periodic bilayer at the energy minimum for case (Λ_AA,Λ_BB) = $Λ_AA, $Λ_BB"

            # Plot the centre of curvature
            scatter!(
                ax2,
                [0.0],
                [0.0],
                color = :red,
                markersize = 10
            )

            Colorbar(fig_bilayer[1,3],colormap = cmap3, colorrange=clims3, label="Couple Stress", width=Relative(0.05),height=Relative(0.6))

            # ax2.title =
            #     "Couple stresses on internal vertices"

            save(datadir(
                "InfiniteBilayer", "coupleStressTests",
                "bilayer_case_(Λ_AA,Λ_BB)=($Λ_AA,$Λ_BB).png"
            ),
            fig_bilayer)
    
            checkHNetwork(R_ext, A, B, computeForces(R_ext,Λ_AA,Λ_BB,Λ_AB,Γ,A,B), Λ_AA,Λ_BB)
    
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
                "Failed $paramInd/$(length(Λ_vec)): ",
                "Λ_AA = $Λ_AA, Λ_BB = $Λ_BB"
            )
            showerror(stderr, err)
            println(stderr)    
    
    
    
    
        end

    end


end



# Vector of Λ values to loop over, [Λ_AA,Λ_BB]
Λ_AB = -0.2
Γ=0.05
nRepeats = 4
Λ_vec = [
        [-0.2,-0.2],
        [-0.3,-0.3],
        [-0.3,-0.25],
        [-0.3,-0.2],
        [-0.3,-0.15],
        [-0.3,-0.1],
        [-0.25,-0.3],
        [-0.25,-0.1],
        [-0.2,-0.339],
        [-0.2,-0.3],
        [-0.15,-0.3],
        [-0.1,-0.3],
        [-0.1,-0.25],
        [-0.339,-0.2],
        ]

isdir(datadir("InfiniteBilayer","coupleStressTests")) ? nothing : mkpath(datadir("InfiniteBilayer","coupleStressTests"))

runBilayer(Λ_vec, Λ_AB, Γ, nRepeats)