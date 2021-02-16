#
#  Visualise.jl
#  VertexModelJL
#
#  Created by Christopher Revell on 16/02/2021.
#
#
#

module Visualise

# Julia packages
using Plots

# Local modules

@inline @views function visualise(Ā,R,nEdges,t,folderName)

   scatter(R[:,1],R[:,2],xlims=(-2,2),ylims=(-2,2),aspect_ratio=:equal,color=:black)

   for i=1:nEdges
      x=findall(x->x>0.5,Ā[i,:])
      plot!([R[x[1],1],R[x[2],1]],[R[x[1],2],R[x[2],2]],color=:black)
   end

   savefig("output/$folderName/plot$t.png")

end

export visualise

end
