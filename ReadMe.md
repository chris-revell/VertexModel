### VertexModel.jl

If running for the first time, open Julia REPL in the directory that contains the VertexModel project directory, not the project directory itself. 
Add the project to local Julia Dev package list: `using Pkg; Pkg.develop(path="./")`. This will allow precompilation to be cached. 

Open a Julia REPL at the project root directory. 

Activate environment with `using Pkg; Pkg.activate(".")`. 

(If running for the first time) Run `Using Pkg; Pkg.instantiate()`

Once `VertexModel.jl` has been precompiled, run function `vertexModel` using default parameters by entering:
`vertexModel()`

Input parameters are optional; use different values other than defaults by running eg. `vertexModel(realTimetMax=12345.0)`

Possible input parameters:
|Variable name:     | Default value:| Explanation:                                                          |
|:------------------|:--------------|:------------------------------------------------------------------    |
|initialSystem      | `"seven"`     | String specifying initial system state                                |
|realTimetMax       | `6.0*86400.0` | Real time maximum system run time /seconds                            |
|realCycleTime      | `86400.0`     | Cell cycle time in seconds                                            |
|γ                  | `0.2`         | Parameter in energy relaxation                                        |
|L₀                 | `0.75`        | Preferred cell perimeter length                                       |
|A₀                 | `1.0`         | Cell preferred area                                                   |
|viscousTimeScale   | `20.0`        | Relaxation rate, approx from Sarah's data.                            |
|pressureExternal   | `0.1`         | External pressure applied isotropically to system boundary            |
|peripheralTension  | `0.0`         | Tension applied to system boundary                                    |
|t1Threshold        | `0.01`        | Edge length at which a T1 transition is triggered                     |
|outputTotal        | `100`         | Number of data output instances per simulation                        |
|outputToggle       | `1`           | Controls whether data are saved from simulation                       |
|frameDataToggle    | `1`           | Controls whether system data are produced at intermediate time points |
|frameImageToggle   | `1`           | Controls whether images are produced at intermediate time points      |
|printToggle        | `1`           | Controls whether updates are printed to command line during runs      |
|videoToggle        | `1`           | Controls whether a movie of the simulation is produced at the end     |
|subFolder          | `""`          | Name of subfolder within data directory in which to store results     |
|solver             | `Tsit5()`     | Name of integration algorithm to use                                  |
|nBlasThreads       | `1`           | Number of threads used by calls to BLAS functions                     |


TODO:
- How do we choose T1 transition edge length threshold? Should this depend on preferred cell perimeter? 
