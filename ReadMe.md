### VertexModel.jl

If running for the first time, open Julia REPL in the directory that contains the VertexModel project directory, not the project directory itself. 
Add the project to local Julia Dev package list: `using Pkg; Pkg.develop(path="./")`. This will allow precompilation to be cached. 

Open a Julia REPL at the project root directory. 

Activate environment with `using Pkg; Pkg.activate(".")`. 

(If running for the first time) Run `Using Pkg; Pkg.instantiate()`

Load package by running `using VertexModel`

Once `VertexModel.jl` has been precompiled, run function `vertexModel` using default parameters by entering:
`vertexModel()`

Input parameters are optional; use different values other than defaults by running eg. `vertexModel(realTimetMax=12345.0)`

Possible input parameters:
|Variable name:     | Default value:| Explanation:                                                               |
|:------------------|:--------------|:---------------------------------------------------------------------------|
|initialSystem      | `"large"`     | String specifying initial system state                                     |
|realTimetMax       | `2.0*86400.0` | Real time maximum system run time /seconds                                 |
|realCycleTime      | `86400.0`     | Cell cycle time in seconds                                                 |
|γ                  | `0.2`         | Parameter in energy relaxation                                             |
|L₀                 | `0.75`        | Preferred cell perimeter length                                            |
|A₀                 | `1.0`         | Cell preferred area                                                        |
|viscousTimeScale   | `20.0`        | Relaxation rate, approx from Sarah's data.                                 |
|pressureExternal   | `0.0`         | External pressure applied isotropically to system boundary                 |
|peripheralTension  | `0.0`         | Tension applied to system boundary                                         |
|t1Threshold        | `0.05`        | Edge length at which a T1 transition is triggered                          |
|solver             | `Tsit5()`     | Name of integration algorithm to use                                       |
|nBlasThreads       | `1`           | Number of threads used by calls to BLAS functions                          |
|subFolder          | `""`          | Name of subfolder within data directory in which to store results          |
|outputTotal        | `100`         | Number of data output instances per simulation                             |
|outputToggle       | `1`           | Controls whether data are saved from simulation                            |
|frameDataToggle    | `1`           | Controls whether system data are produced at intermediate time points      |
|frameImageToggle   | `1`           | Controls whether images are produced at intermediate time points           |
|printToggle        | `1`           | Controls whether updates are printed to command line during runs           |
|videoToggle        | `1`           | Controls whether a movie of the simulation is produced at the end          |
|plotCells          | `1`           | Plot cell polygons in visualisation                                        |
|scatterEdges       | `0`           | Plot scatter points showing cell edge midpoints in visualisation           |
|scatterVertices    | `0`           | Plot scatter points showing vertex locations in visualisation              |
|scatterCells       | `0`           | Plot scatter points showing cell centroids in visualisation                |
|plotForces         | `0`           | Plot arrows showing resultant internal forces on vertices in visualisation |
|setRandomSeed      | `0`           | If non-zero, use this value to set random seed for simulation              |


TODO:
- How do we choose T1 transition edge length threshold? Should this depend on preferred cell perimeter? 
