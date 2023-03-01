### VertexModel.jl

Open Julia REPL at project root directory.

Activate environment with `using Pkg; Pkg.activate(".")` or access package manager with `]` and then enter `activate .` before exiting package manager with backspace. Alternatively, load DrWatson library with `using DrWatson` and enter `@quickactivate`.

If running for the first time, run `Using Pkg; Pkg.instantiate()` or access package manager with `]` and then enter `instantiate` before exiting package manager with backspace.

Load test parameters by entering `includet("scripts/testParameters.jl"). This instantiates variables with values for running the vertex model code. 

These steps can be made easier in future by running `include("loadAll.jl")`.

Once `VertexModel.jl` has been precompiled, run function `vertexModel` using test parameters by entering:
`vertexModel(initialSystem,realTimetMax,realCycleTime,γ,L₀,A₀,viscousTimeScale,dt,pressureExternal,t1Threshold,outputTotal,outputToggle,plotToggle;subFolder="test")`

Input parameters:
initialSystem    (eg. "single")  String specifying initial system state
realTimetMax     (eg. 86400.0 )  Real time maximum system run time /seconds
realCycleTime    (eg. 86400.0 )  Cell cycle time in seconds
γ                (eg. 0.2     )  Parameter in energy relaxation
L₀               (eg. 0.75    )  Preferred cell perimeter length
viscousTimeScale (eg. 20.0    )  Relaxation rate, approx from Sarah's data.
A₀               (eg. 1.0     )  Cell preferred area (1.0 by default)
pressureExternal (eg. 0.2     )  External pressure applied isotropically to system boundary
outputTotal      (eg. 20      )  Number of data outputs
t1Threshold      (eg. 0.01    )  Edge length at which a T1 transition is triggered
outputToggle     (eg. 1       )  Argument controlling whether data are saved from simulation
plotToggle       (eg. 1       )  Argument controlling whether plots are produced from simulation
subFolder        (eg. "Test"  )  Name of subfolder within data directory in which to store results

TODO:
- Control distribution of cell divisions by selecting from random distribution.
- How do we choose T1 transition edge length threshold? Should this depend on preferred cell perimeter? 
