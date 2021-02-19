#
#  Division.jl
#  VertexModelJL
#
#  Created by Christopher Revell on 18/02/2021.
#
#
# Function to update arrays for cell division for all cells that meet division conditions. 

module Division

# Julia packages
using LinearAlgebra

function division!()

    for i=1:nCells
        if CONDITION
            if iseven(cellEdgeCount[i])

            else

            end
        end
    end

end

export division!

end
