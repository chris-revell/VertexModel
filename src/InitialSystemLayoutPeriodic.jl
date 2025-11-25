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
@from "OrderAroundCell.jl" using OrderAroundCell

export initialSystemLayoutPeriodic

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

    len = size(ptsArray,2)
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


function keptEdgeList(ptsArray,triangulation)

    # Function to edit the Delaunay connectivity list such that only edges which dip into the original domain are kept. 
    n = size(ptsArray,2)
    edges = Set{Tuple{Int,Int}}() # Using a set to avoid duplicate edges

    for tri in triangulation.triangles
        # Remap indices to original range 
        tri_mod = [(i-1)%n+1 for i in tri] # mod n and adjusting for 1-based indexing
        # Add edges of the triangle
        push!(edges, Tuple(sort([tri_mod[1], tri_mod[2]])))
        push!(edges, Tuple(sort([tri_mod[2], tri_mod[3]])))
        push!(edges, Tuple(sort([tri_mod[3], tri_mod[1]])))

    end

    # Convert the set of edges to an array 
    edges_array = collect(edges)
    
    return edges_array
end

function keptVerticesList(tessellation,L_x,L_y)
    # Function to only keep vertices which lie within the original periodic domain
    vor_points = SVector{2,Float64}[] # Vector of kept vertices
    kept_indices = Int[]    # Original indices


    for (i,vert) in enumerate(tessellation.polygon_points)
        x, y = vert[1], vert[2]
        if 0.0 <= x <= L_x && 0.0 <= y <= L_y
            push!(vor_points, vert)
            push!(kept_indices, i)
        end
    end
    return kept_indices, vor_points
end



function edges_from_polygons(polygons)
    # Function to build edges from the vertex network
    edges = Set{Tuple{Int,Int}}()

    for poly in polygons
        n = length(poly)
        for i in 1:n-1
            a = poly[i]
            b = poly[i+1]
            push!(edges, (a,b)) # Keep the edge orientation given by polygon order 
        end
    end

    return collect(edges)
end

function buildA(edges, nVerts)
    nEdges = length(edges) # This is okay because edges is a set. 
    A = spzeros(Int, nEdges, nVerts)

    for (ei, (v1,v2)) in enumerate(edges)
        A[ei, v1] = 1
        A[ei, v2] = -1
    end

    return A
end

function buildB(polygons, edges)

    nCells = length(polygons)
    nEdges = length(edges)
    B = spzeros(Int, nCells, nEdges)
    
    # Loop over cells
    for (i,poly) in enumerate(polygons)
        # println(poly)
        # Loop over edges
        n = length(poly)
        for (j,edge) in enumerate(edges)
            v1, v2 = edge 
            # Loop over cell vertices
            if (v1 in poly) && (v2 in poly)
                for k in 1:n-1
                    a = poly[k]
                    b = poly[k+1]
                    # println("a=",a,"b=",b,"v1=",v1,"v2=",v2)
                    if a == v1 && b == v2
                        # edge orientation agrees with cell orientation
                        B[i,j] = 1
                    elseif a == v2 && b == v1
                        B[i,j] = -1
                    end
                end
            end
        end
    end 
    return B
end

function initialSystemLayoutPeriodic(L0_A,L0_B,γ,L_x,L_y)
    # Main function to create periodic initial system layout

    if L0_A == L0_B
        # Compute the roots of the cubic equation in l from the unstressed hexagon area: 
        # Cubic is of the form (9/4)l^3-(sqrt(3)/2 + 6Γ)l + Γ*L0_A. Solve this using the coefficients:
        a,b,c,d = 9/2, 0, (-√(3) + 12*γ), -2*γ*L0_A
        p = Polynomial([d, c, b, a])
        roots_p = roots(p)
        # Only consider real roots
        tol = 1e-10
        real_roots = real.(roots_p[abs.(imag.(roots_p)) .< tol])
        if isempty(real_roots)
            error("No real roots from l cubic")
        else
            l = maximum(real_roots)
        end
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

    # Truncate to N_c of random permutation
    kept = kept[randperm(length(kept))[1:N_c]]

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

    # Only keep vertices within original domain, tracking the old index from tessellation: 
    kept_indices, vor_points = keptVerticesList(tessellation,L_x,L_y)
    println(size(vor_points))
    println("size(kept_indices)",size(kept_indices))

    # Create a dictionary from old vertex indices to new: 
    idx_map = Dict(old => new for (new, old) in enumerate(kept_indices))
    # println(kept_indices)
    # println(enumerate(kept_indices))


    # Convert voronoi_points to an array, R: 
    N_v = length(kept_indices)
    R = zeros(Float64,2,N_v)
    for (i,v) in enumerate(vor_points)
        R[:,i] = v
    end
    R = reinterpret(SVector{2,Float64}, R)

    # The polygon index matches the Delaunay index, so we only keep the first N_c cells which correspond to elements of ptsArray
    N_c = size(ptsArray,2)
    # Generate kept_polygons of same type as tessellation.polygons
    kept_polygons = Vector{Vector{Int}}()
    for i in 1:N_c
        push!(kept_polygons, copy(tessellation.polygons[i]))
    end
    # println("kept_polygons=",kept_polygons)
    for poly in kept_polygons
        
        for (k,old_idx) in enumerate(poly)
            
            if old_idx in keys(idx_map)
                poly[k] = idx_map[old_idx] # map to new index 
            else
                # Find the equivalent point within the domain using mod: 
                vertex = tessellation.polygon_points[old_idx]
                wrapped_vertex = SVector(mod(vertex[1],L_x), mod(vertex[2],L_y))
                new_idx = findfirst(x -> isapprox(x, wrapped_vertex), vor_points)

                if new_idx === nothing
                    println("failed wrapped_vertex=",wrapped_vertex,"failed vertex=",vertex)
                    error("Couldn't find wrapped vertex in vor_points")
                end
                poly[k] = new_idx
            end
        end
        
    end

    # Now determine the cell edges from these kept polygons: 
    # edges = edges_from_polygons(kept_polygons)
    edges = Set{Tuple{Int,Int}}()
    for poly in kept_polygons
        for i in 1:length(poly)-1
            a, b = poly[i], poly[i+1]
            push!(edges, (a,b))  # preserve clockwise order
        end
    end
    edges_array = collect(edges)



    A = buildA(edges_array, N_v)
    B = buildB(kept_polygons, edges_array)

    println("Size A: ", size(A))
    println("Size B: ", size(B))

    
    return A, B, R

end



end
