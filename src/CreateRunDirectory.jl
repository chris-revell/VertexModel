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
    @unpack initialSystem,realTimetMax,γ,λ,A₀,pressureExternal,viscousTimeScale,outputTotal,t1Threshold,realCycleTime,nVerts,nCells,nEdges,L₀,δL,outputInterval,tMax,nonDimCycleTime = params

    # Create directory for run data labelled with current time.
    paramsName = @savename nCells L₀ δL γ t1Threshold realTimetMax
    folderName = "$(paramsName)_$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))"

    # Create frames subfirectory to store system state at each output time
    mkpath(datadir("sims",subFolder,folderName,"frameImages"))
    mkpath(datadir("sims",subFolder,folderName,"frameData"))

    # Store initial system characteristic matrices
    jldsave(datadir("sims",subFolder,folderName,"matricesInitial.jld2");A,B,R)
    writedlm(datadir("sims",subFolder,folderName,"A_initial.csv"), A, ',')
    writedlm(datadir("sims",subFolder,folderName,"B_initial.csv"), B, ',')
    writedlm(datadir("sims",subFolder,folderName,"R_initial.csv"), R, ',')

    return folderName

end

export createRunDirectory

end
