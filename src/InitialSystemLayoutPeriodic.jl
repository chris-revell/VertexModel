#
#  InitialSystemLayoutPeriodic.jl
#  VertexModel
#
#  Created by Charlie Taylor Barca on 11/11/25.
#
#
# Function to create a hexagonal grid of cells. 
# Given number of rows nRows, central row has length nRows, each adjacent row has length nRows-1 etc. 
# Number of cells is then nRows*(nRows-1) - (floor(Int64, nRows/2)+1)*(floor(Int64, nRows/2)+2) + nRows

module InitialSystemLayoutPeriodic

# Julia packages
using LinearAlgebra
using SparseArrays
using StaticArrays
using DrWatson
using FromFile
using DelaunayTriangulation
using FromFile
using Random
using Polynomials 
using Distributions
using Roots
using CairoMakie

@from "SenseCheck.jl" using SenseCheck

function periodic_distance(p, q, Lx, Ly)
    # Function to calculate distance between points p and q, wrapped around the periodic boundaries: 
    dx = abs(p[1] - q[1])
    dy = abs(p[2] - q[2])
    dx = min(dx, Lx - dx)
    dy = min(dy, Ly - dy)
    return sqrt(dx^2 + dy^2)
end

function matern_typeII(parents, marks, rad, Lx, Ly)
    order = sortperm(marks)
    kept = Vector{NTuple{2,Float64}}()

    for idx in order
        p = parents[idx]
        keep = true
        for q in kept
            if periodic_distance(p, q, Lx, Ly) < rad
                keep = false
                break
            end
        end
        keep && push!(kept, p)
    end
    return kept
end

function matern_type2_equation(r, λₚ, λₜ)
    # Equation to solve for the exclusion radius given values of λₚ, λₜ
    return (1 - exp(-λₚ * π * r^2)) / (π * r^2) - λₜ
end

function solve_exclusion_radius(λₚ, λₜ)
    f(r) = matern_type2_equation(r, λₚ, λₜ)

    r_low  = 1e-8              # essentially zero
    r_high = sqrt(1/(π*λₜ))

    r = find_zero(f, (r_low, r_high), Roots.Brent())
    return r
end

function copy_domain_x9(ptsArray,L_x,L_y)
    # Function which copies the cell centres onto a 9x9 grid
    # for Delaunay triangulation on periodic domain 

    len = length(ptsArray[1,:])
    extendedPtsArray = zeros(Float64, (2, 9*len))
    
    extendedPtsArray[:,1:len] = ptsArray 
    extendedPtsArray[:,len+1:2*len] = ptsArray .+ [L_x;0]
    extendedPtsArray[:,2*len+1:3*len] = ptsArray .+ [-L_x;0]
    extendedPtsArray[:,3*len+1:4*len] = ptsArray .+ [0;L_y]
    extendedPtsArray[:,4*len+1:5*len] = ptsArray .+ [0;-L_y]
    extendedPtsArray[:,5*len+1:6*len] = ptsArray .+ [L_x;L_y]
    extendedPtsArray[:,6*len+1:7*len] = ptsArray .+ [-L_x;L_y]
    extendedPtsArray[:,7*len+1:8*len] = ptsArray .+ [L_x;-L_y]
    extendedPtsArray[:,8*len+1:9*len] = ptsArray .+ [-L_x;-L_y]

    return extendedPtsArray
end


function initialSystemLayoutPeriodic(L0_A,L0_B,γ,L_x,L_y)

    

    if L0_A == L0_B
        # Compute the roots of the cubic equation in l from the unstressed hexagon area: 
        # Cubic is of the form (9/4)l^3-(sqrt(3)/2 + 6Γ)l + Γ*L0_A. Solve this using the coefficients:
        a,b,c,d = 9/2, 0, (-sqrt(3) + 12*γ), -2*γ*L0_A
        p = Polynomial([d, c, b, a])
        roots_p = roots(p)
        # Choose the greatest of these: 
        l = maximum(roots_p)
        println("l=",l)
        if l<=0 
            println("Error: negative hexagon sidelength")
        end
        Area_hex = 3*sqrt(3)*l^2/2
        N_c = Int(ceil(L_x*L_y / Area_hex))

        
        # Determine parameters for the Matérn type II process
        λₜ = N_c / (L_x*L_y) # Target intensity 
        λₚ = 10*λₜ # Starting poisson intensity 

        # Solve for exclusion radius: 
        r_ex = solve_exclusion_radius(λₚ, λₜ)

        # Let's try by just making exclusion radius the smallest 
        # distance between hexagon cell centres 
        # r_ex = 2*√3*l/3
        # r_ex = l

        println("r_ex = ",r_ex)
        


        

    else
        error("Cell number calculation for differing L0_A, L0_B not implemented.")
    end


    # Matern type II process to generate periodic cell centres:
    rad = r_ex    
    area = L_x * L_y          # parent intensity guess
    kept = NTuple{2,Float64}[]

    while length(kept) < N_c
        n_parent = rand(Poisson(λₚ * area))
        parents = [(rand()*L_x, rand()*L_y) for _ in 1:n_parent]
        marks   = rand(n_parent)
    
        kept = matern_typeII(parents, marks, rad, L_x, L_y)  
    end

    println("N_c=",N_c)
    println("length(kept)=",length(kept))

    # Rewriting to be in line with InitialSystemLayout.jl
    cellPoints = [SVector(p[1], p[2]) for p in kept]


    xs = [x[1] for x in cellPoints]
    ys = [x[2] for x in cellPoints]

    ptsArray = zeros(Float64, (2, length(cellPoints)))
    for (i, point) in enumerate(cellPoints)
        ptsArray[1, i] = point[1]
        ptsArray[2, i] = point[2]
    end

    extendedPtsArray = copy_domain_x9(ptsArray,L_x,L_y)

    triangulation = triangulate(extendedPtsArray)
    tessellation = voronoi(triangulation, clip=true)


    fig = Figure(resolution=(600,600))
    ax = Axis(fig[1,1], aspect=1, xlabel="x", ylabel="y", title="Voronoi Tessellation")
    
    # Plot Voronoi polygons
    for (gen, vert_indices) in tessellation.polygons
        verts = tessellation.polygon_points[vert_indices]
        xs = [v[1] for v in verts]
        ys = [v[2] for v in verts]
        # Close the polygon
        push!(xs, xs[1])
        push!(ys, ys[1])
        lines!(ax, xs, ys, color=:black)
    end
    
    # Plot original cell centers
    scatter!(ax, extendedPtsArray[1, :], extendedPtsArray[2, :], color=:red, markersize=8)

    # Draw domain boundary
    lines!(ax, [0,L_x,L_x,0,0], [0,0,L_y,L_y,0], color=:blue, linewidth=2)
    
    display(fig)


    # return A, B, R
    return roots_p

end

export initialSystemLayoutPeriodic

end
