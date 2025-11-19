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

function initialSystemLayoutPeriodic(L0_A,L0_B,γ,L_x,L_y)

    

    if L0_A == L0_B
        # Compute the roots of the cubic equation in l from the unstressed hexagon area: 
        # Cubic is of the form (9/4)l^3-(sqrt(3)/2 + 6Γ)l + Γ*L0_A. Solve this using the coefficients:
        a,b,c,d = 9/4, 0, -(sqrt(3)/2 + 6*γ), γ*L0_A
        p = Polynomial([d, c, b, a])
        roots_p = roots(p)
        # Choose the greatest of these: 
        l = maximum(roots_p)
        if l<=0 
            println("Error: negative hexagon sidelength")
        end
        Area_hex = 3*sqrt(3)*l^2/2
        N_c = Int(ceil(L_x*L_y / Area_hex))

        
        

    else
        println("Update cell number calculation for different preferred areas!")
    end


    # Matern type II process to generate periodic cell centres:
    rad = l                      # use hexagon side as hardcore radius
    area = L_x * L_y
    λp = N_c * 3 / area          # parent intensity guess
    kept = NTuple{2,Float64}[]

    while length(kept) < N_c
        n_parent = rand(Poisson(λp * area))
        parents = [(rand()*L_x, rand()*L_y) for _ in 1:n_parent]
        marks   = rand(n_parent)

        kept = matern_typeII(parents, marks, rad, L_x, L_y)

        λp *= 1.5   # increase parent intensity if we failed
    end

    # truncate if too many
    if length(kept) > N_c
        kept = kept[randperm(length(kept))[1:N_c]]
    end

    xs = [x[1] for x in kept]
    ys = [x[2] for x in kept]

    # return A, B, R
    return roots_p

end

export initialSystemLayoutPeriodic

end
