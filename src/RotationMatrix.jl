
# Function to use Rodrigues' Rotation Formula to find the rotation matrix necessary to rotate a 
# vector on the plane of the curved monolayer around a radial vector of the curved monolayer 
# https://mathworld.wolfram.com/RodriguesRotationFormula.html

module RotationMatrix

using LinearAlgebra

# Function takes vector v and returns the rotation matrix in the plane perpendicular to v for a given angle of rotation 
# Direction of rotation follows right-hand rule around given axis vector.
# If no vector passed, defaults to z-aligned unit vector 
# If no angle is given, defaults to 90 degree rotation 
function ϵ(; v=[0.0,0.0,1.0], θ = π/2)
    x, y, z = normalize(v) 
    ω̃ = [
        0 -z y 
        z 0 -x 
        -y x 0
    ]
    
    return (I(3) .+ ω̃.*sin(θ) .+ ω̃*ω̃.*(1-cos(θ)))
end 

function ϵ!(ω̃ ; v=[0.0,0.0,1.0], θ = π/2)
    x, y, z = normalize(v) 
    # ω̃ .= [
    #     0.0 -z y 
    #     z 0.0 -x 
    #     -y x 0.0
    # ]
    
    ω̃[1,1] = 0.0
    ω̃[2,1] = z 
    ω̃[3,1] = -y 
    ω̃[1,2] = -z 
    ω̃[2,2] = 0.0
    ω̃[3,2] = x 
    ω̃[1,3] = y
    ω̃[2,3] = x
    ω̃[3,3] = 0.0
    
    ω̃ .= I(3) .+ ω̃.*sin(θ) .+ ω̃*ω̃.*(1-cos(θ))
    
    return nothing 
end 

export ϵ
export ϵ!

end