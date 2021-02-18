#
#  CreateRunDirectory.jl
#  VertexModelJL
#
#  Created by Christopher Revell on 09/02/2021.
#
#
# Function to create a directory in which to store run results and parameters, with directory name given by currentdate and time

module CreateRunDirectory

# Julia packages
using Dates
using Base.Filesystem
using DelimitedFiles

function createRunDirectory(nCells,nEdges,nVerts,gamma,lamda,tStar,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,A,B,R)

    # Create directory for run data labelled with current time.
    foldername = Dates.format(Dates.now(),"yyyy-mm-dd-HH-MM-SS")
    mkpath("output/$(foldername)")

    # Store system parameters.
    open("output/$(foldername)/conditions.txt","w") do conditionsfile
        println(conditionsfile, "nCells,             $nCells            ")
        println(conditionsfile, "nEdges,             $nEdges            ")
        println(conditionsfile, "nVerts,             $nVerts            ")
        println(conditionsfile, "gamma,              $gamma             ")
        println(conditionsfile, "lamda,              $lamda             ")
        println(conditionsfile, "tStar,              $tStar             ")
        println(conditionsfile, "realTimetMax,       $realTimetMax      ")
        println(conditionsfile, "tMax,               $tMax              ")
        println(conditionsfile, "dt,                 $dt                ")
        println(conditionsfile, "outputInterval,     $outputInterval    ")
        println(conditionsfile, "preferredPerimeter, $preferredPerimeter")
    end

    # Store initial system characteristic matrices
    writedlm("output/$(foldername)/A.txt",A," ")
    writedlm("output/$(foldername)/B.txt",B," ")
    writedlm("output/$(foldername)/R.txt",R," ")

    return foldername

end

export createRunDirectory

end
