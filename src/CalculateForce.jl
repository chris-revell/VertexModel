#
#  CalculateForce.jl
#  VertexModel
#
#  Created by Christopher Revell on 11/02/2021.
#
#
# Function to calculate force vector on vertex k from cell i (Fᵢₖ) for all vertices.

module CalculateForce

# Julia packages
using LinearAlgebra
using StaticArrays
using LoopVectorization
using UnPack
using Plots
using GeometryBasics

function calculateForce!(R,params,matrices)

    @unpack A,B,Ā,B̄,cellTensions,cellPressures,edgeLengths,edgeTangents,F,ϵ = matrices
    @unpack nVerts,nCells,nEdges = params

    fill!(F,@SVector zeros(2))

    # Internal forces
    # NB This iteration could be improved to better leverage sparse arrays
    for k=1:nVerts
        for i=1:nCells
            for j=1:nEdges
                F[k] += 0.5*cellPressures[i]*B[i,j]*Ā[j,k]*(ϵ*edgeTangents[j]) + cellTensions[i]*B̄[i,j]*A[j,k]*edgeTangents[j]/edgeLengths[j]

                # if isnan(F[k][1]) || isnan(F[k][2])
                #     display(R)
                #     plot(xlims=(-0.5,0.5),ylims=(-0.5,0.5),aspect_ratio=:equal,color=:black,markersize=4,markerstroke=:black,dpi=300,size=(1000,1000))
                #     # Scatter vertices
                #     # scatter!(R[:,1],R[:,2])
                #     for i=1:nEdges
                #         x=findall(x->x!=0,Ā[i,:])
                #         # display(x)
                #         plot!([R[x[1]][1],R[x[2]][1]],[R[x[1]][2],R[x[2]][2]],color=:black,linewidth=4)
                #     end
                #     scatter!(Point2f0.(R))
                #     savefig("test.png")
                #     throw()
                # end

                #externalF[k] += boundaryVertices[k]*(0.5*pressureExternal*B[i,j]*Ā[j,k]*(ϵ*edgeTangents[j])) # 0 unless boundaryVertices != 0
            end
        end
    end

end

export calculateForce!

end
