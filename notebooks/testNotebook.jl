
using UnPack
using VertexModel

integ = vertexModel(nRows=5, nCycles=2.0, divisionToggle=0)

(params, matrices) = integ.p 
@unpack A, B = matrices 