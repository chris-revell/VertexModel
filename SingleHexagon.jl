#
#  SingleHexagon.jl
#  VertexModelJL
#
#  Created by Christopher Revell on 15/02/2021.
#
#
#

module SingleHexagon

# Julia packages
using LinearAlgebra

# Local modules

@inline @views function singleHexagon()

    A = [-1.0 1.0  0.0  0.0  0.0  0.0
         0.0  -1.0 1.0  0.0  0.0  0.0
         0.0  0.0  -1.0 1.0  0.0  0.0
         0.0  0.0  0.0  -1.0 1.0  0.0
         0.0  0.0  0.0  0.0  -1.0 1.0
         -1.0 0.0  0.0  0.0  0.0  1.0]

    B = [-1.0 -1.0 -1.0 -1.0 -1.0 1.0]

    R = [0.0 0.0
         0.0 0.0
         0.0 0.0
         0.0 0.0
         0.0 0.0
         0.0 0.0]

    for k=1:6
        R[k,1]= 5.0*(cos((k*π)/3.0))
        R[k,2]= 5.0*(sin((k*π)/3.0))
    end

    return A, B, R

end

export singleHexagon

end
