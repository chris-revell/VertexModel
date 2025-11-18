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

@from "SenseCheck.jl" using SenseCheck

function initialSystemLayoutPeriodic(L0_A,L0_B,γ,L_x,L_y)

    

    if L0_A == L0_B
        # Compute the roots of the cubic equation in l from the unstressed hexagon area: 
        # Cubic is of the form (9/4)l^3-(sqrt(3)/2 + 6Γ)l + Γ*L0_A. Solve this using the coefficients:
        a,b,c,d = 9/4, 0, -(sqrt(3)/2 + 6*γ), γ*L0_A
        p = Polynomial([d, c, b, a])
        roots_p = roots(p)
        # Choose the greatest of these: 
        l = maxiumum(roots_p)
        if l<=0 
            println("Error: negative hexagon sidelength")
        end
        Area_hex = 3*sqrt(3)*l^2/2
        N_c = Int(ceil(L_x*L_y / Area_hex))

        
        

    else
        println("Update cell number calculation for different preferred areas!")
    end

    

    # cellPoints = [SVector(x, 0.0) for x = 1:nRows]
    

    # return A, B, R
    return roots_p

end

export initialSystemLayoutPeriodic

end
