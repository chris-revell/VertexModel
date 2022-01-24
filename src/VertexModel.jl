#
#  VertexModel.jl
#  VertexModel
#
#  Created by Christopher Revell on 31/01/2021.
#
#

module VertexModel

# Julia packages
using LinearAlgebra
using DelimitedFiles
using SparseArrays
using StaticArrays
using Plots
using UnPack
using DrWatson

# Local modules
include("TopologyChange.jl"); using .TopologyChange
include("SpatialData.jl"); using .SpatialData
include("CreateRunDirectory.jl"); using .CreateRunDirectory
include("Visualise.jl"); using .Visualise
include("Initialise.jl"); using .Initialise
include("Iterate.jl"); using .Iterate

# Input parameters:
# initialSystem    (eg. "single")  String specifying initial system state
# realTimetMax     (eg. 86400.0 )  Real time maximum system run time /seconds
# realCycleTime    (eg. 86400.0 )  Cell cycle time in seconds
# γ                (eg. 0.2     )  Parameters in energy relaxation
# λ                (eg. -0.3    )  Parameters in energy relaxation
# viscousTimeScale (eg. 20.0    )  Relaxation rate, approx from Sarah's data.
# dt               (eg. 0.01    )  Non dimensionalised time step
# preferredArea    (eg. 1.0     )  Cell preferred area (1.0 by default)
# pressureExternal (eg. 0.1     )  External pressure applied isotropically to system boundary
# outputTotal      (eg. 20      )  Number of data outputs
# t1Threshold      (eg. 0.01    )  Edge length at which a T1 transition is triggered
# outputToggle     (eg. 1       )  Argument controlling whether data are saved from simulation


function vertexModel(initialSystem,realTimetMax,realCycleTime,γ,λ,viscousTimeScale,dt,preferredArea,pressureExternal,outputTotal,t1Threshold,outputToggle)

    #BLAS.set_num_threads(4)

    # Set up initial system, packaging parameters and matrices for system into params and matrices containers from VertexModelContainers.jl
    params,matrices = initialise(initialSystem,realTimetMax,γ,λ,preferredArea,pressureExternal,dt,viscousTimeScale,outputTotal,t1Threshold,realCycleTime)

    # Extract some variables from containers for use below
    @unpack tMax, outputInterval = params
    @unpack R, tempR, ΔR, cellAges = matrices

    # Initial setup of matrices based on system topology
    topologyChange!(matrices)
    spatialData!(R,params,matrices)

    # Setup output if outputToggle argument == 1
    if outputToggle==1
        # Create fun directory, save parameters, and store directory name for later use.
        folderName = createRunDirectory(params,matrices)
        # Create animation object for visualisation
        anim = Animation()
        visualise(anim,params,matrices)
    end

    t = 1E-8   # Initial time is very small but slightly above 0 to avoid floating point issues with % operator in output interval calculation

    while t<tMax

        # 4 step Runge-Kutta integration
        # 1st step of Runge-Kutta
        iterate!(1,params,matrices)
        # 2nd step of Runge-Kutta
        iterate!(2,params,matrices)
        # 3rd step of Runge-Kutta
        iterate!(2,params,matrices)
        # 4th step of Runge-Kutta
        iterate!(4,params,matrices)

        # Result of Runge-Kutta steps
        R .+= ΔR
        t +=dt
        cellAges .+= dt

        # Visualise system at every output interval
        if t%outputInterval<dt && outputToggle==1
            visualise(anim,params,matrices)
            println("$(t*viscousTimeScale)/$realTimetMax")
        end

    end

    # If outputToggle==1, save animation object as an animated gif and save final system matrices
    if outputToggle==1
        # Store final system characteristic matrices
        jldsave("data/sims/$folderName/matricesFinal.jld2";matrices.A,matrices.B,matrices.R)
        # Save animated gif
        gif(anim, "data/sims/$folderName/animated.gif", fps = 5)
    end

end

export vertexModel

end
