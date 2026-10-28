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

function initialGeometry()
    

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

    u0 = [x10,y10,x20,y20,x30,y30,x50,y50,x60,y60,κ0,λ0]

    return u0
end

function reconstructEquilibrium(x1,y1,x2,y2,x3,y3,x5,y5,x6,y6,κ,λ,A1_vec,A2_vec,d1_vec,d2_vec,d3_vec)
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

    # Vertices around each cell
    cellVertexOrders1 = CircularArray([1,7,8,4,3,2])
    cellVertexOrders2 = CircularArray([3,4,5,6,10,9])

    push!(A1_vec,abs(area(Point{2,Float64}.(R[cellVertexOrders1])))) 
    push!(A2_vec,abs(area(Point{2,Float64}.(R[cellVertexOrders2]))))

    push!(d1_vec,(y1 - y2)/(sqrt(4*x1^2 + (y1-y2)^2)))
    push!(d2_vec,(y5 - y6)/(sqrt(4*x1^2 + (y5-y6)^2)))
    push!(d3_vec,2*y3/(4*x1^2 + 4*y3^2))

    return R, A1_vec, A2_vec, d1_vec, d2_vec, d3_vec
end

function runBilayer()

    ######################################## INITIALISE VECTORS AND PARAMS ########################################

    # parameterLabelVec = ["(O)","(P)","(Q)","(R)","(S)","test","test2"] # Only these because the rest are symmetric

    # Λ_AAVec = [-0.2,-0.1,-0.1,-0.15,-0.3,-0.2,-0.25]
    # Λ_BBVec = [-0.2,-0.25,-0.3,-0.3,-0.3,-0.3,-0.3]
    Λ_AAVec = range(-0.4,0.0,800)
    Λ_BBVec = range(-0.4,0.0,800)

    Λ_AB = -0.2

    # Initialise vectors of each of our minimisation variables 
    κ_grid = zeros(length(Λ_AAVec), length(Λ_BBVec))



    # Initialise values for base case 
    u0 = initialGeometry()

    A1_vec = []
    A2_vec = []
    # Vectors of the ratio of (y1-y2)/l1, (y5-y6)/l8, and 2y_3/l_4, respectively
    d1_vec = []
    d2_vec = []
    d3_vec = []

    for (Aind,Λ_AA) in enumerate(Λ_AAVec)
        for (Bind,Λ_BB) in enumerate(Λ_BBVec)
            


            # Λ_AA = Λ_AAVec[paramInd]
            # Λ_BB = Λ_BBVec[paramInd]

            ######################### INITIALISE BILAYER USING PREVIOUS PARAMS ########################
            
                try 

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

                    u0 = [x10,y10,x20,y20,x30,y30,x50,y50,x60,y60,κ0,λ0]

                    ################################### CALCULATE CELL AREAS ###################################
                    # Reconstruct vertices at equilibrium

                    R, A1_vec, A2_vec, d1_vec, d2_vec, d3_vec = reconstructEquilibrium(x1,y1,x2,y2,x3,y3,x5,y5,x6,y6,κ,λ,A1_vec,A2_vec,d1_vec,d2_vec,d3_vec)
                    
            
                catch err 
                    γ_AB = Λ_AB - (Λ_AA + Λ_BB) / 2

                    # Preserve vector lengths and parameter indexing
                    # push!(x1_vec, NaN)
                    # push!(y1_vec, NaN)
                    # push!(x2_vec, NaN)
                    # push!(y2_vec, NaN)
                    # push!(x3_vec, NaN)
                    # push!(y3_vec, NaN)
                    # push!(x5_vec, NaN)
                    # push!(y5_vec, NaN)
                    # push!(x6_vec, NaN)
                    # push!(y6_vec, NaN)
                    # push!(κ_vec, NaN)
                    # push!(λ_vec, NaN)
                    # push!(γ_AB_vec, γ_AB)

                    push!(A1_vec,NaN)
                    push!(A2_vec,NaN)
                    push!(d1_vec,NaN)
                    push!(d2_vec,NaN)
                    push!(d3_vec,NaN)

                    κ_grid[Aind,Bind] = NaN

            
                    # println(
                    #     stderr,
                    #     "Failed $paramInd/$(length(parameterLabelVec)): ",
                    #     "Λ_AA = $Λ_AA, Λ_BB = $Λ_BB"
                    # )
                    showerror(stderr, err)
                    println(stderr)    
                end


        end
    end

    return κ_grid, Λ_AAVec,Λ_BBVec, A1_vec, A2_vec, d1_vec, d2_vec, d3_vec
end

κ_grid, Λ_AAVec,Λ_BBVec,A1_vec, A2_vec,d1_vec, d2_vec, d3_vec = runBilayer()

A1_grid = reshape(
    A1_vec,
    length(Λ_BBVec),
    length(Λ_AAVec)
)'   # transpose so indices match κ_grid

d1_grid = reshape(
    d1_vec,
    length(Λ_BBVec),
    length(Λ_AAVec)
)'
d2_grid = reshape(
    d2_vec,
    length(Λ_BBVec),
    length(Λ_AAVec)
)'
d3_grid = reshape(
    d3_vec,
    length(Λ_BBVec),
    length(Λ_AAVec)
)'

figA = Figure(size=(700,600))

axA = Axis(
    figA[1,1],
    xlabel=L"\Lambda_{AA}",
    ylabel=L"\Lambda_{BB}",
    title=L"d1"
)

hmA = heatmap!(
    axA,
    Λ_AAVec,
    Λ_BBVec,
    d3_grid,
    colormap=:viridis,
)
# contour!(
#     axA,
#     Λ_AAVec,
#     Λ_BBVec,
#     A1_grid;
#     levels=[0.5],
#     color=:blue,
#     linewidth=3
# )
Colorbar(
    figA[1,2],
    hmA,
    label=L"d_1"
)

# contour!(
#     axA,
#     Λ_AAVec,
#     Λ_BBVec,
#     κ_grid;
#     levels = [0.0],
#     color = :black,
#     linewidth = 3
# )
contour!(
    axA,
    Λ_AAVec,
    Λ_BBVec,
    d1_grid;
    levels = [0.5],
    color = :pink,
    linewidth = 3
)
contour!(
    axA,
    Λ_AAVec,
    Λ_BBVec,
    d1_grid;
    levels = [0.5],
    color = :green,
    linewidth = 3
)
contour!(axA,
    Λ_AAVec,
    Λ_BBVec,
    d1_grid;
    levels = [0.761513637058961],
    color = :red,
    linewidth = 3)

display(figA)

zeroPairs = findNonDiagonalZeroContour(
    κ_grid,
    Λ_AAVec,
    Λ_BBVec;
    diagonalTolerance = 2e-3
)


###################### PLOTTING INFORMATION ALONG κ=0 ######################
figure = Figure(size=(700,600))

axis = Axis(
    figure[1,1],
    xlabel=L"Λ_{AA}",
    ylabel=L"d3" 
)

for (pairInd,zeroPair) in enumerate(zeroPairs)
    Λ_AA = zeroPair[1]
    Λ_BB = zeroPair[2]

    A1_vec = []
    A2_vec = []
    # Vectors of the ratio of (y1-y2)/l1, (y5-y6)/l8, and 2y_3/l_4, respectively
    d1_vec = []
    d2_vec = []
    d3_vec = []

    u0 = initialGeometry()

    f(u) = energyCalculation(u, Λ_AA, Λ_BB, Λ_AB)

    try 
        result = optimize(
            f,
            u0,
            LBFGS();
            autodiff = :forward,
        )

        u = Optim.minimizer(result)

        x1,y1,x2,y2,x3,y3,x5,y5,x6,y6,κ,λ = u

        R, A1_vec, A2_vec, d1_vec, d2_vec, d3_vec = reconstructEquilibrium(x1,y1,x2,y2,x3,y3,x5,y5,x6,y6,κ,λ,A1_vec,A2_vec,d1_vec,d2_vec,d3_vec)

        scatter!(axis, [Λ_AA], d3_vec, color =:red, markersize = 5)
    catch err 

    end

    # Now plot the diagonal equivalent 
    Λ_BB = Λ_AA 

    A1_vec = []
    A2_vec = []
    # Vectors of the ratio of (y1-y2)/l1, (y5-y6)/l8, and 2y_3/l_4, respectively
    d1_vec = []
    d2_vec = []
    d3_vec = []

    u0 = initialGeometry()

    try 
        result = optimize(
            f,
            u0,
            LBFGS();
            autodiff = :forward,
        )

        u = Optim.minimizer(result)

        x1,y1,x2,y2,x3,y3,x5,y5,x6,y6,κ,λ = u

        R, A1_vec, A2_vec, d1_vec, d2_vec, d3_vec = reconstructEquilibrium(x1,y1,x2,y2,x3,y3,x5,y5,x6,y6,κ,λ,A1_vec,A2_vec,d1_vec,d2_vec,d3_vec)

        scatter!(axis, [Λ_AA], d3_vec, color =:black, markersize = 5)
    catch err 

    end


end
display(figure)



# Find the difference between the values to see where it crosses the diagonal;
d = zeros(length(zeroPairs))
d .= (zeroPair[1]-zeroPair[2] for zeroPair in zeroPairs)

ΛAA_cross = []
ΛBB_cross = []
for i in 1:length(d)-1
    if d[i] * d[i+1] < 0

        t = -d[i] / (d[i+1] - d[i])

        ΛAA_cross = zeroPairs[i][1] + t*(zeroPairs[i+1][1]-zeroPairs[i][1])
        ΛBB_cross = zeroPairs[i][2] + t*(zeroPairs[i+1][2]-zeroPairs[i][2])

        println("Crossing at:")
        println("(Λ_AA,Λ_BB)=($ΛAA_cross,$ΛBB_cross)")
    end
end

# # Look to see where the nonlinear contour crosses the line γ_AB = 0
# Λ_AB = -0.2
# γ_AB_line = zeros(length(zeroPairs))
# γ_AB_line .= ((Λ_AB - (zeroPair[1]+zeroPair[2])/2) for zeroPair in zeroPairs)

# for i in 1:length(γ_AB_line)-1
#     if γ_AB_line[i] * γ_AB_line[i+1] < 0

#         t = -γ_AB_line[i] / (γ_AB_line[i+1] - γ_AB_line[i])

#         ΛAA_cross_γ = zeroPairs[i][1] + t*(zeroPairs[i+1][1]-zeroPairs[i][1])
#         ΛBB_cross_γ = zeroPairs[i][2] + t*(zeroPairs[i+1][2]-zeroPairs[i][2])

#         println("Crossing at:")
#         println("(Λ_AA,Λ_BB)=($ΛAA_cross_γ,$ΛBB_cross_γ)")
#     end
# end


# Find the minimum distance: 
findall(x -> x == minimum(d), d)

# println("zeroAAs = ", [zeroPair[1] for zeroPair in zeroPairs])
# println("zeroBBs = ", [zeroPair[2] for zeroPair in zeroPairs])

fig1 = Figure(size=(700,600))

ax1 = Axis(
    fig1[1,1],
    xlabel=L"Λ_{AA}",
    ylabel=L"Λ_{BB}"
)

scatter!(
    ax1,
    [zeroPair[1] for zeroPair in zeroPairs],
    [zeroPair[2] for zeroPair in zeroPairs],
    color = :black,
    markersize = 10,
    label = "κ=0"
)

display(fig1)

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
    colormap=:RdBu
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

# Define the value at the floppy boundary: 
Γ = 0.05
Λ_star = -2*2*sqrt(2)*3^(1/4)*Γ 
println(Λ_star)
lines!(
    ax,
    [Λ_star,Λ_star],
    [-0.4,0.0],
    color =:green,
    linestyle =:dash,
    linewidth = 2
)
lines!(
    ax,
    [-0.4,0.0],
    [Λ_star,Λ_star],
    color =:green,
    linestyle =:dash,
    linewidth = 2
)


for (crossInd, Λ_AA) in enumerate(ΛAA_cross)
    Λ_BB = ΛBB_cross[crossInd]
    Λ_AB = -0.2

    # Λ_BB = -0.1
    # Λ_AA = -0.1


    u0 = initialGeometry()

    f(u) = energyCalculation(u, Λ_AA, Λ_BB, Λ_AB)

    result = optimize(
        f,
        u0,
        LBFGS();
        autodiff = :forward,
    )

    u = Optim.minimizer(result)

    x1,y1,x2,y2,x3,y3,x5,y5,x6,y6,κ,λ = u

    d3 = 2*y3/(4*x1^2 + 4*y3^2)

    # R, A1, A2, d1, d2, d3 = reconstructEquilibrium(x1,y1,x2,y2,x3,y3,x5,y5,x6,y6,κ,λ,A1_vec,A2_vec,d1_vec,d2_vec,d3_vec)
    println("d3_vec = ",d3)
end



scatter!(ax, [ΛAA_cross], [ΛBB_cross],color =:cyan, markersize=10)

Colorbar(fig[1,2], hm, label = L"κ", width = 20)

display(fig)
save(datadir("InfiniteBilayer","bilayer_κ_heatmap.png"), fig)


