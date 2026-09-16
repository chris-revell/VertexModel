# This notebook will run the vertex model and unpack the incidence matrices for use in the REPL 

using UnPack
using VertexModel

# Run an initial simulation and output the integrator object to the integ variable
integ = vertexModel(nRows=5, 
    nCycles=2.0, 
    divisionToggle=false,
)

# Unpack properties of the simulation from the integ variable
(params, matrices) = integ.p 
@unpack A, B = matrices 
R = reinterpret(SVector{2,Float64}, integ.u)

# Here you can perform operations to change the state of the system 

# Now run another simulation starting from the R, A, B state
integ2 = vertexModel(initialSystem="argument", 
    R_in = R, 
    A_in = A, 
    B_in = B,nCycles=2.0, 
    divisionToggle=false,
)