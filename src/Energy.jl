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

    @unpack cellEnergies,
        cellAreas,
        cellPerimeters = matrices
    @unpack nCells,
        A₀,
        L₀,
        γ = params

    energyTotal = sum(cellAreas.*(log.(cellAreas).-1)+(γ*L₀).*cellPerimeters.*(log.(cellPerimeters./L₀).-1))

    return energyTotal

end

export energy

end
