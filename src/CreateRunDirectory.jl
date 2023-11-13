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
using DrWatson

function createRunDirectory(R,params,matrices,subFolder)

    @unpack A,B = matrices
    @unpack initialSystem,realTimetMax,γ,λ,A₀,pressureExternal,viscousTimeScale,outputTotal,t1Threshold,realCycleTime,nVerts,nCells,nEdges,L₀,outputInterval,tMax,nonDimCycleTime = params

    # Create directory for run data labelled with current time.
    paramsName = @savename L₀ γ t1Threshold realTimetMax
    folderName = "$(paramsName)_$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))"
    # Create frames subfirectory to store system state at each output time
    mkpath(datadir("sims",subFolder,folderName,"frameImages"))
    mkpath(datadir("sims",subFolder,folderName,"frameData"))

    # Store system parameters.
    # open(datadir("sims",subFolder,folderName,"conditions.txt"),"w") do conditionsFile
    #     println(conditionsFile, "initialSystem,     $initialSystem")
    #     println(conditionsFile, "realTimetMax,      $realTimetMax")
    #     println(conditionsFile, "realCycleTime,     $realCycleTime")
    #     println(conditionsFile, "γ,                 $γ")
    #     println(conditionsFile, "λ,                 $λ")
    #     println(conditionsFile, "viscousTimeScale,  $viscousTimeScale")
    #     println(conditionsFile, "A₀,                $A₀")
    #     println(conditionsFile, "pressureExternal,  $pressureExternal")
    #     println(conditionsFile, "outputTotal,       $outputTotal")
    #     println(conditionsFile, "t1Threshold,       $t1Threshold")
    #     println(conditionsFile, "nVerts,            $nVerts")
    #     println(conditionsFile, "nCells,            $nCells")
    #     println(conditionsFile, "nEdges,            $nEdges")
    #     println(conditionsFile, "L₀,                $L₀")
    #     println(conditionsFile, "outputInterval,    $outputInterval")
    #     println(conditionsFile, "tMax,              $tMax")
    #     println(conditionsFile, "nonDimCycleTime,   $nonDimCycleTime")
    # end

    return folderName

end

export createRunDirectory

end
