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
using LoopVectorization
using Plots
using UnPack
using DrWatson
@quickactivate

# Local modules
include("TopologyChange.jl"); using .TopologyChange
include("CreateRunDirectory.jl"); using .CreateRunDirectory
include("Visualise.jl"); using .Visualise
include("Initialise.jl"); using .Initialise
include("Iterate.jl"); using .Iterate


function vertexModel(initialSystem,realTimetMax,γ,λ,tStar,dt,preferredArea,pressureExternal,outputTotal,t1Threshold,outputToggle)

    #BLAS.set_num_threads(4)

    # Input parameters
    # initialSystem    (eg. "single")  String specifying initial system state
    # realTimetMax     (eg. 10000.0 )  Real time maximum system run time /seconds
    # γ                (eg. 0.2     )  Parameters in energy relaxation
    # λ                (eg. -0.3    )  Parameters in energy relaxation
    # tStar            (eg. 20.0    )  Relaxation rate, approx from Sarah's data.
    # dt               (eg. 0.01    )  Non dimensionalised time step
    # preferredArea    (eg. 1.0     )  Cell preferred area (1.0 by default)
    # pressureExternal (eg. 0.1     )  External pressure applied isotropically to system boundary
    # outputTotal      (eg. 20      )  Number of data outputs
    # t1Threshold      (eg. 0.01    )  Edge length at which a T1 transition is triggered
    # outputToggle     (eg. 1       )  Argument controlling whether data are saved from simulation

    cellCycleTime = 24

    # Set up initial system, packaging parameters and matrices for system into params and matrices containers from VertexModelContainers.jl
    params,matrices = initialise(initialSystem,realTimetMax,γ,λ,preferredArea,pressureExternal,dt,tStar,outputTotal,t1Threshold,cellCycleTime)

    # Initial setup of matrices based on system topology
    topologyChange!(matrices)

    # Extract some variables from containers
    @unpack tMax, outputInterval = params
    @unpack R, tempR, ΔR, cellAges = matrices

    # Setup output if outputToggle argument == 1
    if outputToggle==1
        # Create fun directory, save parameters, and store directory name for later use.
        folderName = createRunDirectory(params,matrices)
        # Create animation object for visualisation
        anim = Animation()
        # Run first visualisation step to save initial system state
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
            println("$t/$tMax")
        end

    end

    # If outputToggle==1, save animation object as an animated gif
    outputToggle==1 ? gif(anim, "data/sims/$folderName/animated.gif", fps = 10) : nothing

end

export vertexModel

end
