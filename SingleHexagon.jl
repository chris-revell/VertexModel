#
#  SingleHexagon.jl
#  VertexModelJL
#
#  Created by Christopher Revell on 15/02/2021.
#
#
# Function to create a starting system of one hexagonal cell.

module SingleHexagon

# Julia packages
using LinearAlgebra

# Local modules

@inline @views function singleHexagon()

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

    return A, B, R

end

export singleHexagon

end
