#
#  CreateRunDirectory.jl
#  VertexModel
#
#  Created by Christopher Revell on 09/02/2021.
#
#
# Function to create a directory in which to store simulations results and parameters, with directory name given by current date and time

module CreateRunDirectory

# Julia packages
using Dates
using Base.Filesystem
using DelimitedFiles
using UnPack

function createRunDirectory(R,nVerts, nCells, nEdges, γ, λ, preferredPerimeter, preferredArea, pressureExternal, dt, outputInterval, tStar, realTimetMax, tMax,A,B)

    # @unpack A,B = matrices
    # @unpack nVerts, nCells, nEdges, γ, λ, preferredPerimeter, preferredArea, pressureExternal, dt, outputInterval, tStar, realTimetMax, tMax = params

    # Create directory for run data labelled with current time.
    foldername = Dates.format(Dates.now(),"yyyy-mm-dd-HH-MM-SS")
    mkpath("data/sims/$(foldername)")

    # Store system parameters.
    open("data/sims/$(foldername)/conditions.txt","w") do conditionsfile
        println(conditionsfile, "nCells,             $nCells            ")
        println(conditionsfile, "nEdges,             $nEdges            ")
        println(conditionsfile, "nVerts,             $nVerts            ")
        println(conditionsfile, "γ,                  $γ                 ")
        println(conditionsfile, "λ,                  $λ                 ")
        println(conditionsfile, "tStar,              $tStar             ")
        println(conditionsfile, "realTimetMax,       $realTimetMax      ")
        println(conditionsfile, "tMax,               $tMax              ")
        println(conditionsfile, "dt,                 $dt                ")
        println(conditionsfile, "outputInterval,     $outputInterval    ")
        println(conditionsfile, "preferredPerimeter, $preferredPerimeter")
        println(conditionsfile, "preferredArea,      $preferredArea     ")
    end

    # Store initial system characteristic matrices
    writedlm("data/sims/$(foldername)/A.txt",A," ")
    writedlm("data/sims/$(foldername)/B.txt",B," ")
    writedlm("data/sims/$(foldername)/R.txt",R," ")

    return foldername

end

export createRunDirectory

end
