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

function createRunDirectory(params,matrices)

    @unpack A,B,R = matrices
    @unpack initialSystem,realTimetMax,γ,λ,preferredArea,pressureExternal,dt,viscousTimeScale,outputTotal,t1Threshold,realCycleTime = params

    # Create directory for run data labelled with current time.
    folderName = Dates.format(Dates.now(),"yyyy-mm-dd-HH-MM-SS")
    mkpath("data/sims/$(folderName)")

    # Store system parameters.
    open("data/sims/$(folderName)/conditions.txt","w") do conditionsfile
        println(conditionsfile, "initialSystem,     $initialSystem"
        println(conditionsfile, "realTimetMax,      $realTimetMax"
        println(conditionsfile, "γ,                 $γ"
        println(conditionsfile, "λ,                 $λ"
        println(conditionsfile, "preferredArea,     $preferredArea"
        println(conditionsfile, "pressureExternal,  $pressureExternal"
        println(conditionsfile, "dt,                $dt"
        println(conditionsfile, "viscousTimeScale,  $viscousTimeScale"
        println(conditionsfile, "outputTotal,       $outputTotal"
        println(conditionsfile, "t1Threshold,       $t1Threshold"
        println(conditionsfile, "realCycleTime,     $realCycleTime"
    end

    # Store initial system characteristic matrices
    writedlm("data/sims/$(folderName)/Ainitial.txt",A," ")
    writedlm("data/sims/$(folderName)/Binitial.txt",B," ")
    writedlm("data/sims/$(folderName)/Rinitial.txt",R," ")

    return folderName

end

export createRunDirectory

end
