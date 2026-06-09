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

function initialEdgeLength(γ,L₀)
    a,b,c,d = 9/2, 0, (-√(3) + 12*γ), -2*γ*L₀
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
    return l
end

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




function initialSystemLayoutPeriodic(γ,L_x,L_y,Λ_AA,Λ_BB,Λ_AB,Area_A_ratio)
    # Main function to create periodic initial system layout

        # Compute effective preferred perimeters for isolated A- or B-cells
        L0_A = -Λ_AA/(2*γ)
        L0_B = -Λ_BB/(2*γ)
        println("L0_A=", L0_A)
        println("L0_B=", L0_B)

        # Compute initial edge lengths for t1 thresholds: 
        l_AA = initialEdgeLength(γ, -Λ_AA/(2*γ))
        l_BB = initialEdgeLength(γ, -Λ_BB/(2*γ))
        l_AB = initialEdgeLength(γ, -Λ_AB/(2*γ))

        # Desired ratio of Area_A : Area_B
        # Area_A_ratio = 0.3
        Area_B_ratio = 1.0 - Area_A_ratio

        Area_A = Area_A_ratio * (L_x * L_y)
        Area_B = Area_B_ratio * (L_x * L_y)

        # Compute the roots of the cubic equation in l from the unstressed hexagon area: 
        # Cubic is of the form (9/4)l^3-(sqrt(3)/2 + 6Γ)l + Γ*L0_A. Solve this using the coefficients:

        
        Area_hex_A = 3*sqrt(3)*l_AA^2/2
        N_cA = Int(ceil(Area_A/ Area_hex_A))

        Area_hex_B = 3*sqrt(3)*l_BB^2/2
        N_cB = Int(ceil(Area_B/ Area_hex_B))
        
        # Determine parameters for the Matérn type II process
        λₜA = N_cA / (L_x*L_y*Area_A_ratio) # Target intensity 
        λₚA = 20*λₜA # Starting poisson intensity 

        λₜB = N_cB / (L_x*L_y*(Area_B_ratio)) # Target intensity 
        λₚB = 20*λₜB # Starting poisson intensity


        # Solve for exclusion radius: 
        r_exA = solve_exclusion_radius(λₚA, λₜA)
        r_exB = solve_exclusion_radius(λₚB, λₜB)

        # Matern type II process to generate periodic cell centres:
        rad = r_exA +r_exB   
        area = L_x * L_y          # parent intensity guess
        kept = NTuple{2,Float64}[]

        while length(kept) < N_cA
            n_parent = rand(Poisson(λₚA * area))
            parents = [(rand()*L_x, rand()*L_y) for _ in 1:n_parent]
            marks   = rand(n_parent)
        
            kept = matern_typeII(parents, marks, rad, L_x, L_y)  
        end

        # Truncate to N_c of random permutation
        kept = kept[randperm(length(kept))[1:N_cA]]

        # Now add B cells:
        # rad = r_exB
        while length(kept) < N_cA + N_cB
            n_parent = rand(Poisson(λₚB * area))
            parents = [(rand()*L_x, rand()*L_y) for _ in 1:n_parent]
            marks   = rand(n_parent)
        
            new_kept = matern_typeII(parents, marks, rad, L_x, L_y)  
            append!(kept, new_kept)
        end

        # Truncate to N_c of random permutation
        kept = kept[randperm(length(kept))[1:N_cA+N_cB]]

    
        # Rewriting to be in line with InitialSystemLayout.jl
        cellPoints = [SVector(p[1], p[2]) for p in kept]

        ε = min(r_exA/5, 0.1)
        for i in 1:length(cellPoints)
            cellPoints[i] += SVector(randn()*ε, randn()*ε)
            cellPoints[i] = SVector(mod(cellPoints[i][1], L_x),
                                mod(cellPoints[i][2], L_y))
        end



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

        # Create a dictionary from old vertex indices to new: 
        idx_map = Dict(old => new for (new, old) in enumerate(kept_indices))


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
        orderedPairs = Set{Tuple{Int,Int}}()
        for poly in kept_polygons
            for i in 1:length(poly)-1
                a, b = poly[i], poly[i+1]
                push!(orderedPairs, (min(a,b),max(a,b)))  # preserve clockwise order
            end
        end

        nVerts = length(R)
        nEdges = length(orderedPairs)
        nCells = length(cellPoints)

        # Construct A matrix mapping tessellation edges to tessellation vertices 
        A = spzeros(Int64, nEdges, nVerts)
        for (edgeIndex, vertices) in enumerate(orderedPairs)
            A[edgeIndex, vertices[1]] = 1
            A[edgeIndex, vertices[2]] = -1
        end

        B = spzeros(Int64, nCells, nEdges)
        for c = 1:nCells
            for i = 2:length(kept_polygons[c])
                vertexLeading = kept_polygons[c][i-1]  # Leading with respect to *clockwise* direction around cell
                vertexTrailing = kept_polygons[c][i]
                # Find index of edge connecting these vertices 
                edge = (findall(x -> x != 0, @view A[:, vertexLeading])∩findall(x -> x != 0, @view A[:, vertexTrailing]))[1]
                if A[edge, vertexLeading] > 0
                    B[c, edge] = 1
                else
                    B[c, edge] = -1
                end
            end
        end


        
        return A, B, R, N_cA, l_AA, l_BB, l_AB

    

    

end



end
