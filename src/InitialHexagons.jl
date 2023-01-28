#
#  InitialHexagons.jl
#  VertexModel
#
#  Created by Christopher Revell on 15/02/2021.
#
#
# Function to create an initial system of 1 or 3 hexagonal cells.

module InitialHexagons

# Julia packages
using LinearAlgebra
using SparseArrays
using StaticArrays

indexLoop(a,N) = (N+a-1)%(N)+1

function initialHexagons(n)

    # Argument of n=1 produces a single initial hexagonal cell
    if n=="one"

        ATmp= [-1.0 1.0 0.0 0.0 0.0 0.0
             0.0 -1.0 1.0 0.0 0.0 0.0
             0.0 0.0 -1.0 1.0 0.0 0.0
             0.0 0.0 0.0 -1.0 1.0 0.0
             0.0 0.0 0.0 0.0 -1.0 1.0
             1.0 0.0 0.0 0.0 0.0 -1.0]
        A = sparse(ATmp)

        B = -1.0.*ones(1,6)

        R = Array{SVector{2,Float64}}(undef,6)
        for k=1:6
            R[k] = SVector{2}([cos((k*π)/3.0),sin((k*π)/3.0)])
        end

        # Set initial cell areas to 1.0
        R .*= 1.0/(2*sin(π/3.0)*(1+cos(π/3.0)))

    elseif n=="three"

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

        R = Array{SVector{2,Float64}}(undef,13)
        a = sin(π/3)
        b = cos(π/3)
        R[1]  = SVector{2}([-1.0,-2*a])
        R[2]  = SVector{2}([0.0,-2*a])
        R[3]  = SVector{2}([b,-a])
        R[4]  = SVector{2}([0.0,0.0])
        R[5]  = SVector{2}([-1.0-b,-a])
        R[6]  = SVector{2}([-1.0,0.0])
        R[7]  = SVector{2}([-1.0-b,a])
        R[8]  = SVector{2}([-1.0,2.0*a])
        R[9]  = SVector{2}([0.0,2.0*a])
        R[10] = SVector{2}([b,a])
        R[11] = SVector{2}([1.0+b,-a])
        R[12] = SVector{2}([1.0+2*b,0.0])
        R[13] = SVector{2}([1.0+b,a])
        # Set initial cell areas to 1.0
        R .*= 1.0/(2*sin(π/3.0)*(1+cos(π/3.0)))

    elseif n=="seven"
        A = spzeros(Int64,18,12)
        B = spzeros(Int64,7, 18)
        R = Array{SVector{2,Float64}}(undef,12)

        B[1,1]  = 1
        B[1,2]  = 1
        B[1,3]  = 1
        B[1,4]  = 1
        B[1,5]  = 1
        B[1,6]  = 1
        B[2,1]  = -1
        B[2,7]  = 1
        B[2,8]  = 1
        B[2,9]  = 1
        B[3,2]  = -1
        B[3,9]  = -1
        B[3,10] = 1
        B[3,11] = 1
        B[4,3]  = -1
        B[4,11] = -1
        B[4,12] = 1
        B[4,13] = 1
        B[5,4]  = -1
        B[5,13] = -1
        B[5,14] = 1
        B[5,15] = 1
        B[6,5]  = -1
        B[6,15] = -1
        B[6,16] = 1
        B[6,17] = 1
        B[7,6]  = -1
        B[7,7]  = -1
        B[7,17] = -1
        B[7,18] = 1

        A[1,1] = -1
        A[1,2] = 1
        A[2,2] = -1
        A[2,3] = 1
        A[3,3] = -1
        A[3,4] = 1
        A[4,4] = -1
        A[4,5] = 1
        A[5,5] = -1
        A[5,6] = 1
        A[6,6] = -1
        A[6,1] = 1
        A[7,1] = -1
        A[7,7] = 1
        A[8,7] = -1
        A[8,8] = 1
        A[9,8] = -1
        A[9,2] = 1
        A[10,8] = -1
        A[10,9] = 1
        A[11,9] = -1
        A[11,3] = 1
        A[12,9] = -1
        A[12,10] = 1
        A[13,10] = -1
        A[13,4] = 1
        A[14,10] = -1
        A[14,11] = 1
        A[15,11] = -1
        A[15,5] = 1
        A[16,11] = -1
        A[16,12] = 1
        A[17,12] = -1
        A[17,6] = 1
        A[18,12] = -1
        A[18,7] = 1

        for k=1:6
            R[k] = SVector{2}([cos((-k*π)/3.0),sin((-k*π)/3.0)])
            R[k+6] = 3.0.*SVector{2}([cos((-k*π)/3.0),sin((-k*π)/3.0)])
        end

        # # Set initial cell areas to 1.0
        R .*= 1.0/(2*sin(π/3.0)*(1+cos(π/3.0)))



    end

    return A, B, R

end

export initialHexagons

end
