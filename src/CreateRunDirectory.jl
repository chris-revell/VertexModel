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
using JLD2 

function createRunDirectory(params,matrices)

    @unpack A,B,R = matrices
    @unpack initialSystem,realTimetMax,γ,λ,preferredArea,pressureExternal,dt,viscousTimeScale,outputTotal,t1Threshold,realCycleTime,nVerts,nCells,nEdges,preferredPerimeter,outputInterval,tMax,nonDimCycleTime = params

    # Create directory for run data labelled with current time.
    folderName = Dates.format(Dates.now(),"yyyy-mm-dd-HH-MM-SS")
    mkpath("data/sims/$(folderName)")

    # Store system parameters.
    open("data/sims/$(folderName)/conditions.txt","w") do conditionsFile
        println(conditionsFile, "initialSystem,     $initialSystem")
        println(conditionsFile, "realTimetMax,      $realTimetMax")
        println(conditionsFile, "realCycleTime,     $realCycleTime")
        println(conditionsFile, "γ,                 $γ")
        println(conditionsFile, "λ,                 $λ")
        println(conditionsFile, "viscousTimeScale,  $viscousTimeScale")
        println(conditionsFile, "dt,                $dt")
        println(conditionsFile, "preferredArea,     $preferredArea")
        println(conditionsFile, "pressureExternal,  $pressureExternal")
        println(conditionsFile, "outputTotal,       $outputTotal")
        println(conditionsFile, "t1Threshold,       $t1Threshold")
        println(conditionsFile, "nVerts,            $nVerts")
        println(conditionsFile, "nCells,            $nCells")
        println(conditionsFile, "nEdges,            $nEdges")
        println(conditionsFile, "preferredPerimeter,$preferredPerimeter")
        println(conditionsFile, "outputInterval,    $outputInterval")
        println(conditionsFile, "tMax,              $tMax")
        println(conditionsFile, "nonDimCycleTime,   $nonDimCycleTime")
    end

    jldsave("data/sims/$(folderName)/params.jld2";params)

    # Store initial system characteristic matrices
    jldsave("data/sims/$(folderName)/matricesInitial.jld2";A,B,R)

    return folderName

end

export createRunDirectory

end
