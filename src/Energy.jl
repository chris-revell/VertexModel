#
#  Energy.jl
#  VertexModel
#
#  Created by Christopher Revell on 31/01/2022.
#
# Function to calculate system energy

module Energy

# Julia packages
using LinearAlgebra
using UnPack


function energy(params,matrices)

    @unpack cellAreas,
        cellA₀s,
        cellPerimeters,
        cellL₀s,
        μ,
        Γ = matrices
    @unpack energyModel = params
    
    
    # Quadratic energy
    if energyModel == "quadratic"
        energyTotal = sum(μ.*(0.5 .* (cellAreas .- cellA₀s).^2 .+ 0.5 .* Γ .* (cellPerimeters .- cellL₀s).^2))
    end
    
    

    return energyTotal
end

export energy


end
