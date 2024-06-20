#
#  Energy.jl
#  VertexModel
#
#  Created by Christopher Revell on 31/01/2022.
#
#
# Function to calculate system energy

module Energy

# Julia packages
using LinearAlgebra
using UnPack

function energy(params,matrices)

    @unpack cellAreas,
        cellPerimeters = matrices
    @unpack nCells,
        A₀,
        L₀,
        γ = params

    energyTotal = sum(0.5.*((cellAreas.-A₀).^2+γ.*(cellPerimeters.-L₀).^2))

    return energyTotal

end

export energy

end
