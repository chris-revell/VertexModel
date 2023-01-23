#
#  SenseCheck.jl
#  VertexModel
#
#  Created by Christopher Revell on 18/01/2023.
#
#
# Function to calculate check for nonzero values in B*A

module SenseCheck

# Julia packages
using LinearAlgebra
using SparseArrays

function senseCheck(A, B; marker="")

    test = B*A
    dropzeros!(test)
    if length(findnz(test)[1]) > 0
        throw("Non-zero values in BA: $(marker)")
    else
        return
    end

end

export senseCheck

end
