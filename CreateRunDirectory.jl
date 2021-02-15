#
#  CreateRunDirectory.jl
#  VertexModelJL
#
#  Created by Christopher Revell on 09/02/2021.
#
#
# Function to create a directory in which to store run results, with directory name given by current time

module CreateRunDirectory

using Dates
using Base.Filesystem
using DelimitedFiles

function createRunDirectory(nCells,nEdges,nVerts,A,B,R)

    # Create directory for run data labelled with current time.
    foldername = Dates.format(Dates.now(),"yyyy-mm-dd-HH-MM-SS")
    mkpath("output/$(foldername)")

    # Store system parameters.
    open("output/$(foldername)/conditions.txt","w") do conditionsfile
        println(conditionsfile,"nCells,      $nCells")
        println(conditionsfile,"nEdges,      $nEdges")
        println(conditionsfile,"nVerts,      $nVerts")
    end

    writedlm("output/$(foldername)/A.txt",A," ")
    writedlm("output/$(foldername)/B.txt",B," ")
    writedlm("output/$(foldername)/R.txt",R," ")

    return foldername
end

export createRunDirectory

end
