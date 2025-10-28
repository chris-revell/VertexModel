### VertexModel.jl

Open a Julia REPL at the project root directory. 

Activate environment with `using Pkg; Pkg.activate(".")`, or `] activate .`. 

(If running for the first time) Run `Using Pkg; Pkg.instantiate()` or `] instantiate`.

Load package by running `using VertexModel`

Once `VertexModel.jl` has been precompiled, run function `vertexModel` using default parameters by entering:
`vertexModel()`

Input parameters are optional; use different values other than defaults by running eg. `vertexModel(nRows=7)`

Possible input parameters:

|Variable name:           | Default value:| Explanation:                                                                           |
|:------------------------|:-----------------------|:---------------------------------------------------------------------------------------|
| `initialSystem`         | `"new"`                | Specify whether initial state is generated afresh or imported from file or argument. Possible values `"new"`, `"argument"`, or path to data file. 
| `nRows`                 | `9`                    | Size of initial system state if generating afresh
| `nCycles`               | `1`                    | Number of cell cycles in simulation run time
| `realCycleTime`         | `86400.0`              | Dimensional cell cycle time
| `realTimetMax`          | `nCycles*realCycleTime`| Dimensional simulation run time 
| `γ`                     | `0.2`                  | Parameter in energy relaxation
| `L₀`                    | `0.75`                 | Preferred cell perimeter length
| `A₀`                    | `1.0`                  | Preferred cell area, almost always =1.0
| `viscousTimeScale`      | `1000.0` | Relaxation rate
| `pressureExternal`      | `0.0` | Pressure perpendicular to system boundary applied to boundary vertices
| `peripheralTension`     | `0.0` | Tension parallel to system boundary applied to boundary vertices
| `t1Threshold`           | `0.05` | Edge length below which a T1 transition occurs
| `divisionToggle`        | `1` | Toggle controlling whether cells divide in the simulation. 
| `solver`                | `Tsit5()` | ODE solver algorithm; requires loading DifferentialEquations.jl or OrdinaryDiffEq.jl in REPL to pass solver as argument
| `nBlasThreads`          | `1` | Number of BLAS parallel threads
| `subFolder`             | `""` | Subfolder into which data are saved 
| `outputTotal`           | `100` | Number of data outputs per run
| `outputToggle`          | `1` | Controls whether anything is output from the simulation at all
| `frameDataToggle`       | `1` | Controls whether `.jld2` data files are saved for each output
| `frameImageToggle`      | `1` | Controls whether images are saved for each output (often not helpful on HPC)
| `printToggle`           | `1` | Controls whether updates are printed to the REPL during a run (often not helpful on HPC)
| `videoToggle`           | `1` | Controls whether a video of the simulation is saved  (often not helpful on HPC)
| `plotCells`             | `1` | Controls whether cell polygons are plotted in image output
| `scatterEdges`          | `0` | Controls whether edge midpoints are plotted in image output
| `scatterVertices`       | `0` | Controls whether vertex positions are plotted in image output
| `scatterCells`          | `0` | Controls whether cell midpoints are plotted in image output
| `plotForces`            | `0` | Controls whether force vectors are plotted in image output
| `plotEdgeMidpointLinks` | `0` | Controls whether lines connecting adjacent edge midpoints are plotted in image output
| `randomSeed`            | `0` | Used to specify a random seed for reproducibility; generated from current time if `0` is passed.
| `abstol`                | `1e-7,` | Absolute tolerance for ODE solver (see DifferentialEquations.jl)
| `reltol`                | `1e-4` | Relative tolerance for ODE solver (see DifferentialEquations.jl)
| `energyModel`           | `"log"` | Energy model used; options are `"log"` or `"quadratic"`
| `vertexWeighting`       | `1` | Controls whether drag proportional to vertex area is applied to each vertex
| `R_in`                  | `spzeros(2)` | Vertex position matrix passed at command line; only used if `initialSystem=="argument"`
| `A_in`                  | `spzeros(2)` | Incidence matrix passed at command line; only used if `initialSystem=="argument"`
| `B_in`                  | `spzeros(2)` | Incidence matrix passed at command line; only used if `initialSystem=="argument"`
