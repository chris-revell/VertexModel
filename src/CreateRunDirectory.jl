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
using LibGit2

function createRunDirectory(params,subFolder)

    @unpack initialSystem,
        realTimetMax,
        γ,
        λ,
        A₀,
        pressureExternal,
        viscousTimeScale,
        outputTotal,
        t1Threshold,
        realCycleTime,
        nVerts,
        nCells,
        nEdges,
        L₀,
        outputInterval,
        tMax,
        nonDimCycleTime = params

    repo = LibGit2.GitRepo(projectdir())
    branchname = LibGit2.shortname(LibGit2.head(repo))
    # Create directory for run data labelled with current time.
    paramsName = @savename nCells L₀ γ realCycleTime realTimetMax
    # folderName = joinpath(branchname, subFolder, "$(branchname)_$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
    folderName = joinpath("sims", branchname, subFolder, "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)")
    # Create frames subdirectory to store system state at each output time
    mkpath(datadir(folderName, "frameImages"))
    mkpath(datadir(folderName, "frameData"))
    divName= @savename L₀ γ realCycleTime 
    open(io -> writedlm(io, ["time" "cellIndex" "daughterIndex_1" "daughterIndex_2" "cellPosition_x" "cellPosition_y" "cellLineage" "parentGeneration" "longAxis_x" "longaxis_y"], ','), datadir(folderName,"Division_$(divName).csv"), "w") # write
    return folderName

end

export createRunDirectory

end
