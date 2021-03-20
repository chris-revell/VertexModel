#
#  InitialHexagons.jl
#  VertexModelJL
#
#  Created by Christopher Revell on 15/02/2021.
#
#
# Function to create a starting system hexagonal cells.

module InitialHexagons

# Julia packages
using LinearAlgebra
using SparseArrays

# Local modules

@inline @views function initialHexagons(n)

    if n==1
        A = [-1.0 1.0 0.0 0.0 0.0 0.0
             0.0 -1.0 1.0 0.0 0.0 0.0
             0.0 0.0 -1.0 1.0 0.0 0.0
             0.0 0.0 0.0 -1.0 1.0 0.0
             0.0 0.0 0.0 0.0 -1.0 1.0
             1.0 0.0 0.0 0.0 0.0 -1.0]

        B = -1.0.*ones(1,6)

        R = zeros(6,2)

        for k=1:6
            R[k,1] = cos((k*π)/3.0)
            R[k,2] = sin((k*π)/3.0)
        end
        # Set initial cell areas to 1.0
        R.*=1.0/(2*sin(π/3.0)*(1+cos(π/3.0)))

    elseif n==3

        A = spzeros(Int64,15,13)
        A[1,1] = -1
        A[1,2] = 1
        A[2,2] = -1
        A[2,3] = 1
        A[3,3] = -1
        A[3,4] = 1
        A[4,5] = -1
        A[4,1] = 1
        A[5,6] = -1
        A[5,5] = 1
        A[6,7] = -1
        A[6,6] = 1
        A[7,8] = -1
        A[7,7] = 1
        A[8,9] = -1
        A[8,8] = 1
        A[9,10] = -1
        A[9,9] = 1
        A[10,4] = -1
        A[10,10] = 1
        A[11,13] = -1
        A[11,10] = 1
        A[12,12] = -1
        A[12,13] = 1
        A[13,11] = -1
        A[13,12] = 1
        A[14,3] = -1
        A[14,11] = 1
        A[15,4] = -1
        A[15,6] = 1

        B = spzeros(Int64,3,15)
        B[1,1] = -1
        B[1,2] = -1
        B[1,3] = -1
        B[1,4] = -1
        B[1,5] = -1
        B[1,15] = -1
        B[2,3] = 1
        B[2,14] = -1
        B[2,13] = -1
        B[2,12] = -1
        B[2,11] = -1
        B[2,10] = 1
        B[3,15] = 1
        B[3,10] = -1
        B[3,9] = -1
        B[3,8] = -1
        B[3,7] = -1
        B[3,6] = -1

        R = zeros(13,2)
        a = sin(π/3)
        b = cos(π/3)
        R[1,:] .= [-1.0,-2*a]
        R[2,:] .= [0.0,-2*a]
        R[3,:] .= [b,-a]
        R[4,:] .= [0.0,0.0]
        R[5,:] .= [-1.0-b,-a]
        R[6,:] .= [-1.0,0.0]
        R[7,:] .= [-1.0-b,a]
        R[8,:] .= [-1.0,2.0*a]
        R[9,:] .= [0.0,2.0*a]
        R[10,:] .= [b,a]
        R[11,:] .= [1.0+b,-a]
        R[12,:] .= [1.0+2*b,0.0]
        R[13,:] .= [1.0+b,a]
        # Set initial cell areas to 1.0
        R.*=1.0/(2*sin(π/3.0)*(1+cos(π/3.0)))

    end

    return A, B, R

end

export initialHexagons

end
