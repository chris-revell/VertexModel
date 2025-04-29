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

# Energy per Cowley et al. 2024 Section 2a
𝒰(θ) = θ*(log(θ)-1.0)
Uᵢ(Aᵢ, A₀, Lᵢ, L₀, μᵢ, Γᵢ) = μᵢ*(𝒰(Aᵢ/A₀) + Γᵢ*L₀^2*𝒰(Lᵢ/L₀))

function energy_log(params,matrices)

    @unpack cellAreas,
        cellA₀s,
        cellPerimeters,
        cellL₀s,
        μ,
        Γ = matrices
    
    energyTotal = 0.0
    for i = 1:nCells
        # energyTotal += 0.5 * (cellAreas[i] - A₀)^2 + 0.5 * γ * (cellPerimeters[i] - L₀)^2
        energyTotal += Uᵢ.(cellAreas[i], cellA₀s[i], cellPerimeters[i], cellL₀s[i], μ[i], Γ[i])
    end
    return energyTotal
end


function energy_quadratic(params,matrices)

    @unpack cellAreas,
        cellA₀s,
        cellPerimeters,
        cellL₀s,
        μ,
        Γ = matrices
    
    energyTotal= sum(μ.*(0.5 .* (cellAreas .- cellA₀s).^2 .+ 0.5 .* Γ .* (cellPerimeters .- cellL₀s).^2))

    return energyTotal
end

export energy_log, energy_quadratic

end
