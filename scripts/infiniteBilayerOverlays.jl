
# Code to overlay similar parameter cases

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

function runBilayer(parameterLabelVec,Λ_AAVec,Λ_BBVec,Λ_AB)

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

    fig_bilayer = Figure(size=(900, 700))
    ax1 = Axis(
        fig_bilayer[1, 1],
        aspect = DataAspect()
    )

    hidedecorations!(ax1)
    hidespines!(ax1)

    # Number of periodic copies on either side of the fundamental unit
    nRepeats = 0
    text!(
        ax1,
        0.0,
        0.0,
        text = "O",
        offset = (6, 6),
        align = (:left, :bottom),
        color = :red
    )

    ax1.title =
        "Overlay of periodic bilayer geometries for different parameter cases"

    colVec = ["red","black"]


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

                push!(γ_AB_vec, Λ_AB - (Λ_AA+Λ_BB)/2)

                ########################## Plot the bilayer at equilibrium ##########################

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

                # Vertices around each cell
                cellVertexOrders1 = CircularArray([1,7,8,4,3,2])
                cellVertexOrders2 = CircularArray([3,4,5,6,10,9])

                polygon1 = makeBilayerPolygon(R,cellVertexOrders1)
                polygon2 = makeBilayerPolygon(R,cellVertexOrders2)

                if plotMultipleImages
                    

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
                            colourA = RGB(1,1,1)
                            colourB = RGB(1,1,1)
                            lineWidth = 0.5
                        else
                            colourA = RGBA(1,1,1)
                            colourB = RGBA(1,1,1)
                            lineWidth = 0.5
                        end

                        poly!(
                            ax1,
                            polygon1_image,
                            color = RGBA(1,1,1,0),
                            strokecolor = colVec[paramInd],
                            strokewidth = lineWidth
                        )

                        poly!(
                            ax1,
                            polygon2_image,
                            color = RGBA(1,1,1,0),
                            strokecolor = colVec[paramInd],
                            strokewidth = lineWidth
                        )
                    end

                    # Plot the centre of curvature
                    scatter!(
                        ax1,
                        [0.0],
                        [0.0],
                        color =colVec[paramInd],
                        markersize = 10
                    )

                    

                    
                else
                    

                    fig_bilayer = Figure(resolution=(600,600))
                    ax1 = Axis(fig_bilayer[1,1],aspect=1)
                    hidedecorations!(ax1)
                    hidespines!(ax1)

                    poly!(ax1, polygon1, strokecolor=(:black, 1.0), strokewidth=1)
                    poly!(ax1, polygon2, strokecolor=(:black, 1.0), strokewidth=1)

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
                    save(datadir("InfiniteBilayer","multipleRuns","bilayer_case_$(parameterLabelVec[paramInd])_geometry.png"), fig_bilayer)

                end

                save(datadir(
                            "InfiniteBilayer", "overlayPlots",
                            "bilayer_case_$(parameterLabelVec[paramInd])_periodic_geometry.png"
                        ),
                        fig_bilayer
                    )
        
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


# Standard parameter sets
 
Λ_AAVec = [-0.15,-0.13]
Λ_BBVec = [-0.15,-0.15]

Λ_AB = -0.2

# parameterLabelVec = ["(O)","(P)","(Q)","(R)","(S)","test","test2"]
parameterLabelVec = []
for (ind,AAval) in enumerate(Λ_AAVec)
    BBval = Λ_BBVec[ind]
    push!(parameterLabelVec, ["(Λ_{AA},Λ_{BB}) = ($AAval,$BBval)"])
end


x1_vec, y1_vec, x2_vec, y2_vec, x3_vec, y3_vec, x5_vec, y5_vec, x6_vec, y6_vec, κ_vec, λ_vec,γ_AB_vec = runBilayer(parameterLabelVec,Λ_AAVec,Λ_BBVec,Λ_AB)


