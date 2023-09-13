### VertexModel.jl

Open Julia REPL at project root directory.

(When running for the first time) Add project to local Julia Dev package list: `using Pkg; Pkg.develop(path="./")`. This will allow precompilation to be cached. 

Activate environment with `using Pkg; Pkg.activate(".")`. 

(If running for the first time) Run `Using Pkg; Pkg.instantiate()`

Once `VertexModel.jl` has been precompiled, run function `vertexModel` using default parameters by entering:
`vertexModel()`

Input parameters are optional; use different values other than defaults by running eg. `vertexModel(realTimetMax=12345.0)`

Possible input parameters:
|Variable name:     | Default value:| Explanation:                                                      |
|:------------------|:--------------|:------------------------------------------------------------------|
|initialSystem      | `"seven"`     | String specifying initial system state                            |
|realTimetMax       | `6.0*86400.0` | Real time maximum system run time /seconds                        |
|realCycleTime      | `86400.0`     | Cell cycle time in seconds                                        |
|γ                  | `0.2`         | Parameter in energy relaxation                                    |
|L₀                 | `0.75`        | Preferred cell perimeter length                                   |
|A₀                 | `1.0`         | Cell preferred area                                               |
|viscousTimeScale   | `20.0`        | Relaxation rate, approx from Sarah's data.                        |
|pressureExternal   | `0.1`         | External pressure applied isotropically to system boundary        |
|peripheralTension  | `0.0`         | Tension applied to system boundary                                |
|t1Threshold        | `0.01`        | Edge length at which a T1 transition is triggered                 |
|outputTotal        | `100`         | Number of data output instances per simulation                    |
|outputToggle       | `1`           | Controls whether data are saved from simulation                   |
|plotToggle         | `1`           | Controls whether plots are produced from simulation               |
|subFolder          | `""`          | Name of subfolder within data directory in which to store results |
|solver             | `Tsit5()`     | Name of integration algorithm to use                              |

TODO:
- Control distribution of cell divisions by selecting from random distribution.
- How do we choose T1 transition edge length threshold? Should this depend on preferred cell perimeter? 
