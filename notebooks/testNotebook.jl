# This notebook will run the vertex model and unpack the incidence matrices for use in the REPL 

using UnPack
using VertexModel

integ = vertexModel(nRows=5, nCycles=2.0, divisionToggle=0)

(params, matrices) = integ.p 
@unpack A, B = matrices 
R = reinterpret(SVector{2,Float64}, integ.u)