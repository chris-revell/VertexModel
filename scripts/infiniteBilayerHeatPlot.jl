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

    # Vector of edge tensions for each edge: 
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

function findNonDiagonalZeroContour(
    κ_grid,
    Λ_AAVec,
    Λ_BBVec;
    diagonalTolerance = 1e-4,
    zeroTolerance = 1e-10
)
    zeroPairs = Tuple{Float64,Float64}[]

    # For each fixed Λ_AA, search along the Λ_BB direction
    for Aind in eachindex(Λ_AAVec)

        Λ_AA = Λ_AAVec[Aind]

        for Bind in 1:(length(Λ_BBVec) - 1)

            κ1 = κ_grid[Aind, Bind]
            κ2 = κ_grid[Aind, Bind + 1]

            Λ_BB1 = Λ_BBVec[Bind]
            Λ_BB2 = Λ_BBVec[Bind + 1]

            # Ignore failed or non-finite optimisation results
            if !isfinite(κ1) || !isfinite(κ2)
                continue
            end

            Λ_BB_zero = nothing

            # Grid point is already approximately zero
            if abs(κ1) < zeroTolerance
                Λ_BB_zero = Float64(Λ_BB1)

            # κ changes sign between neighbouring grid points
            elseif κ1 * κ2 < 0
                # Linear interpolation of the κ = 0 location
                Λ_BB_zero =
                    Float64(Λ_BB1 - κ1 * (Λ_BB2 - Λ_BB1) / (κ2 - κ1))
            end

            if Λ_BB_zero !== nothing
                # Remove the diagonal solution Λ_AA = Λ_BB
                if abs(Λ_AA - Λ_BB_zero) > diagonalTolerance
                    push!(
                        zeroPairs,
                        (Float64(Λ_AA), Λ_BB_zero)
                    )
                end
            end
        end
    end

    return zeroPairs
end

function runBilayer()

    ######################################## INITIALISE VECTORS AND PARAMS ########################################

    # parameterLabelVec = ["(O)","(P)","(Q)","(R)","(S)","test","test2"] # Only these because the rest are symmetric

    # Λ_AAVec = [-0.2,-0.1,-0.1,-0.15,-0.3,-0.2,-0.25]
    # Λ_BBVec = [-0.2,-0.25,-0.3,-0.3,-0.3,-0.3,-0.3]
    Λ_AAVec = range(-0.3,-0.1,500)
    Λ_BBVec = range(-0.3,-0.1,500)

    # Initialise vectors of each of our minimisation variables 
    κ_grid = zeros(length(Λ_AAVec), length(Λ_BBVec))



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

    for (Aind,Λ_AA) in enumerate(Λ_AAVec)
        for (Bind,Λ_BB) in enumerate(Λ_BBVec)
            


            # Λ_AA = Λ_AAVec[paramInd]
            # Λ_BB = Λ_BBVec[paramInd]

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

            κ_grid[Aind,Bind] = κ

            # Set initial values for next run: 
            λ0 = λ
            κ0 = κ
            x10 = x1
            y10 = y1
            x20 = x2
            y20 = y2
            x30 = x3
            y30 = y3
            x50 = x5
            y50 = y5
            x60 = x6
            y60 = y6


        end
    end

    return κ_grid, Λ_AAVec,Λ_BBVec
end

κ_grid, Λ_AAVec,Λ_BBVec = runBilayer()

zeroPairs = findNonDiagonalZeroContour(
    κ_grid,
    Λ_AAVec,
    Λ_BBVec;
    diagonalTolerance = 2e-3
)

println("zeroAAs = ", [zeroPair[1] for zeroPair in zeroPairs])

println("zeroBBs = ", [zeroPair[2] for zeroPair in zeroPairs])
fig = Figure(size=(700,600))

ax = Axis(
    fig[1,1],
    xlabel=L"Λ_{AA}",
    ylabel=L"Λ_{BB}",
    title=L"κ"
)

hm = heatmap!(
    ax,
    Λ_AAVec,
    Λ_BBVec,
    κ_grid,
    colormap=:RdBu,
)

contour!(
    ax,
    Λ_AAVec,
    Λ_BBVec,
    κ_grid;
    levels = [0.0],
    color = :black,
    linewidth = 3
)

Colorbar(fig[1,2], hm, label = L"κ", width = 20)

display(fig)
save(datadir("InfiniteBilayer","bilayer_κ_heatmap.png"), fig)
