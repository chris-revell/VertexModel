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
    @unpack nCells,A₀,L₀,γ = params

    energyTotal = 0.0
    for i=1:nCells
        energyTotal += 0.5*(cellAreas[i]-A₀)^2 + 0.5*γ*(cellPerimeters[i]-L₀)^2
    end

    return energyTotal

end

export energy

end
