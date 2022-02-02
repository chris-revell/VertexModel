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

    @unpack cellEnergies,cellAreas,cellPerimeters = matrices
    @unpack nCells,preferredArea,preferredPerimeter,γ = params

    energyTotal = 0.0
    for i=1:nCells
        energyTotal += 0.5*(cellAreas[i]-preferredArea)^2 + 0.5*γ*(cellPerimeters[i]-preferredPerimeter)^2
    end
    # cellEnergies .= 0.5.*(cellAreas.-preferredArea).^2 .+ 0.5*γ.*(cellPerimeters.-preferredPerimeter).^2
    # energyTotal = sum(cellEnergies)

    return energyTotal

end

export energy

end
