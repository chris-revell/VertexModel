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
    if n==1

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

    elseif n==7
        A = spzeros(Int64,30,24)
        B = spzeros(Int64,7, 30)
        R = Array{SVector{2,Float64}}(undef,24)

        B[1,1] = 1
        B[1,2] = 1
        B[1,3] = 1
        B[1,4] = 1
        B[1,5] = 1
        B[1,6] = 1
        B[2,2] = -1
        B[2,19] = 1
        B[2,20] = 1
        B[2,21] = 1
        B[2,22] = 1
        B[2,23] = 1
        B[3,3] = -1
        B[3,15] = 1
        B[3,16] = 1
        B[3,17] = 1
        B[3,18] = 1
        B[3,19] = -1
        B[4,4] = -1
        B[4,11] = 1
        B[4,12] = 1
        B[4,13] = 1
        B[4,14] = 1
        B[4,15] = -1
        B[5,5] = -1
        B[5,7] = 1
        B[5,8] = 1
        B[5,9] = 1
        B[5,10] = 1
        B[5,11] = -1
        B[6,6] = -1
        B[6,7] = -1
        B[6,27] = 1
        B[6,28] = 1
        B[6,29] = 1
        B[6,30] = 1
        B[7,1] = -1
        B[7,27] = -1
        B[7,26] = 1
        B[7,25] = 1
        B[7,24] = 1
        B[7,23] = -1

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
        A[7,7] = -1
        A[7,6] = 1
        A[8,8] = -1
        A[8,7] = 1
        A[9,9] = -1
        A[9,8] = 1
        A[10,10] = -1
        A[10,9] = 1
        A[11,10] = -1
        A[11,5] = 1
        A[12,11] = -1
        A[12,10] = 1
        A[13,12] = -1
        A[13,11] = 1
        A[14,13] = -1
        A[14,12] = 1
        A[15,13] = -1
        A[15,4] = 1
        A[16,14] = -1
        A[16,13] = 1
        A[17,15] = -1
        A[17,14] = 1
        A[18,16] = -1
        A[18,15] = 1
        A[19,16] = -1
        A[19,3] = 1
        A[20,17] = -1
        A[20,16] = 1
        A[21,18] = -1
        A[21,17] = 1
        A[22,19] = -1
        A[22,18] = 1
        A[23,2] = -1
        A[23,19] = 1
        A[24,20] = -1
        A[24,19] = 1
        A[25,21] = -1
        A[25,20] = 1
        A[26,22] = -1
        A[26,21] = 1
        A[27,1] = 1
        A[27,22] = -1
        A[28,23] = -1
        A[28,22] = 1
        A[29,24] = -1
        A[29,23] = 1
        A[30,7] = -1
        A[30,24] = 1

        for k=1:6
            R[k] = SVector{2}([cos((-k*π)/3.0),sin((-k*π)/3.0)])
        end

        edgeLength1 = 1.0
        edgeLength2 = 2.0*sin(π/3.0)

        R[22] = R[1] .+ SVector{2}([cos((-1*π)/3.0),sin((-1*π)/3.0)]).*edgeLength1
        R[21] = R[1] .+ SVector{2}([cos((-1.5*π)/3.0),sin((-1.5*π)/3.0)]).*edgeLength2
        R[20] = R[2] .+ SVector{2}([cos((-1.5*π)/3.0),sin((-1.5*π)/3.0)]).*edgeLength2

        R[19] = R[2] .+ SVector{2}([cos((-2*π)/3.0),sin((-2*π)/3.0)]).*edgeLength1
        R[18] = R[2] .+ SVector{2}([cos((-2.5*π)/3.0),sin((-2.5*π)/3.0)]).*edgeLength2
        R[17] = R[3] .+ SVector{2}([cos((-2.5*π)/3.0),sin((-2.5*π)/3.0)]).*edgeLength2

        R[16] = R[3] .+ SVector{2}([cos((-3*π)/3.0),sin((-3*π)/3.0)]).*edgeLength1
        R[15] = R[3] .+ SVector{2}([cos((-3.5*π)/3.0),sin((-3.5*π)/3.0)]).*edgeLength2
        R[14] = R[4] .+ SVector{2}([cos((-3.5*π)/3.0),sin((-3.5*π)/3.0)]).*edgeLength2

        R[13] = R[4] .+ SVector{2}([cos((-4*π)/3.0),sin((-4*π)/3.0)]).*edgeLength1
        R[12] = R[4] .+ SVector{2}([cos((-4.5*π)/3.0),sin((-4.5*π)/3.0)]).*edgeLength2
        R[11] = R[5] .+ SVector{2}([cos((-4.5*π)/3.0),sin((-4.5*π)/3.0)]).*edgeLength2

        R[10] = R[5] .+ SVector{2}([cos((-5*π)/3.0),sin((-5*π)/3.0)]).*edgeLength1
        R[9] = R[5] .+ SVector{2}([cos((-5.5*π)/3.0),sin((-5.5*π)/3.0)]).*edgeLength2
        R[8] = R[6] .+ SVector{2}([cos((-5.5*π)/3.0),sin((-5.5*π)/3.0)]).*edgeLength2

        R[7] = R[6] .+ SVector{2}([cos((-6*π)/3.0),sin((-6*π)/3.0)]).*edgeLength1
        R[24] = R[6] .+ SVector{2}([cos((-6.5*π)/3.0),sin((-6.5*π)/3.0)]).*edgeLength2
        R[23] = R[1] .+ SVector{2}([cos((-6.5*π)/3.0),sin((-6.5*π)/3.0)]).*edgeLength2




        # Set initial cell areas to 1.0
        R .*= 1.0/(2*sin(π/3.0)*(1+cos(π/3.0)))



    end

    return A, B, R

end

export initialHexagons

end
